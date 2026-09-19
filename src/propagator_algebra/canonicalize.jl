"""
    advanced_to_retarded(x, prefactor)

Apply the transformation ``G^A(y,y)=-G^R(y,y)`` to equal-position contractions.
"""
function advanced_to_retarded(
    x::Vector{Contraction{S}}, prefactor::Number
) where {S<:Statistics}
    ff(c::Contraction{S}) = is_advanced(c) && same_position(c)
    adv_idx = findall(ff, x)
    if isempty(adv_idx)
        return x, prefactor
    end
    x′ = copy(x)
    for i in adv_idx
        prefactor *= -1
        x′[i] = adjoint(x[i])
    end
    return x′, prefactor
end

function sort_by_position_and_type(p::Contraction)::Float64
    if is_out(p)
        return -Inf
    elseif is_in(p)
        return Inf
    else
        i, j = integer_positions(p)
        type = Int(propagator_type(p...))
        return float(pairing(i, j) * 4 + type)
    end
end
function sort_by_position_and_type(
    p::Tuple{Field{S},Field{S}}
)::Float64 where {S<:Statistics}
    return sort_by_position_and_type(Contraction(p))
end
sort_by_position_and_type(p::Edge)::Float64 =
    sort_by_position_and_type(Contraction(fields(p)))

@inline function field_state_color(f::Field)::UInt64
    derivative = UInt64(derivative_multiindex(f).orders)
    orientation_bits = UInt64(Int(orientation(f))) << 63
    keldysh_bits = UInt64(Int(keldysh_index(f))) << 62
    regularisation_bits = UInt64(Int(regularisation(f)) + 1) << 60
    return orientation_bits | keldysh_bits | regularisation_bits | derivative
end

@inline function field_color(f::Field)
    return (name(f), slots(field_indices(f)), field_state_color(f))
end

@inline function propagator_color(c::Contraction)
    return (field_color(c.out), field_color(c.in), Int(propagator_type(c...)))
end
@inline function propagator_color(e::Edge)
    return (field_color(e.out), field_color(e.in), Int(propagator_type(e)))
end

function canonicalization_positions(vs)
    result = Position[]
    sizehint!(result, 2 * length(vs))
    for item in vs
        for p in positions(item)
            p in result || push!(result, p)
        end
    end
    sort!(result)
    return result
end

# Uniform-color detection is on the canonicalization hot path. Cache the reference field
# families, packed state colors, and propagator type once, then compare each remaining edge
# against those concrete values. Derivative decoration remains part of the packed state while
# the derivative-free path avoids rebuilding nested color tuples or decoding the reference
# fields on every comparison.
function uniform_coloring(vs)
    isempty(vs) && return true
    reference = first(vs)
    reference_out_family = field_family(reference.out)
    reference_in_family = field_family(reference.in)
    reference_out_state = field_state_color(reference.out)
    reference_in_state = field_state_color(reference.in)
    reference_type = propagator_type(reference)

    @inbounds for i in (firstindex(vs) + 1):lastindex(vs)
        item = vs[i]
        isequal(field_family(item.out), reference_out_family) || return false
        isequal(field_family(item.in), reference_in_family) || return false
        field_state_color(item.out) == reference_out_state || return false
        field_state_color(item.in) == reference_in_state || return false
        propagator_type(item) === reference_type || return false
    end
    return true
end

function simple_position_pairs(vs)
    for i in eachindex(vs)
        item_positions = positions(vs[i])
        for j in firstindex(vs):(i - 1)
            isequal(positions(vs[j]), item_positions) && return false
        end
    end
    return true
end

function position_labels(graph_positions::Vector{Position})
    labels = Vector{Int}(undef, length(graph_positions))
    for (i, p) in enumerate(graph_positions)
        labels[i] = if is_out(p)
            1
        elseif is_in(p)
            2
        else
            3
        end
    end
    return labels
end

@inline position_vertex(graph_positions::Vector{Position}, p::Position) =
    searchsortedfirst(graph_positions, p)

# Build the direct position graph and report whether every directed position pair
# was unique. Detecting multiplicity while inserting edges avoids an O(E^2)
# duplicate-edge pre-scan on the canonicalization hot path.
function _make_simple_NautyDiGraph(vs, graph_positions::Vector{Position})
    graph = NautyGraphs.NautyDiGraph(
        length(graph_positions); vertex_labels=position_labels(graph_positions)
    )
    simple = true
    for item in vs
        out, in = positions(item)
        source = position_vertex(graph_positions, out)
        target = position_vertex(graph_positions, in)
        simple &= !Graphs.has_edge(graph, source, target)
        Graphs.add_edge!(graph, source, target)
    end
    return graph, simple
end

function make_simple_NautyDiGraph(vs, graph_positions::Vector{Position})
    graph, _ = _make_simple_NautyDiGraph(vs, graph_positions)
    return graph
end

function propagator_colors(vs)
    C = typeof(propagator_color(first(vs)))
    colors = C[]
    sizehint!(colors, length(vs))
    for item in vs
        color = propagator_color(item)
        color in colors || push!(colors, color)
    end
    sort!(colors)
    return colors
end

function make_colored_NautyDiGraph(vs, graph_positions::Vector{Position})
    colors = propagator_colors(vs)
    npositions = length(graph_positions)
    labels = Vector{Int}(undef, npositions + length(vs))
    copyto!(labels, 1, position_labels(graph_positions), 1, npositions)

    for (i, item) in enumerate(vs)
        color_index = searchsortedfirst(colors, propagator_color(item))
        labels[npositions + i] = 3 + color_index
    end

    graph = NautyGraphs.NautyDiGraph(length(labels); vertex_labels=labels)
    for (i, item) in enumerate(vs)
        out, in = positions(item)
        edge_vertex = npositions + i
        Graphs.add_edge!(graph, position_vertex(graph_positions, out), edge_vertex)
        Graphs.add_edge!(graph, edge_vertex, position_vertex(graph_positions, in))
    end
    return graph
end

const GC_PHYSICAL_MAX_VERTICES = 64

"""Reusable released-GC scratch for physical propagator canonicalization."""
mutable struct PhysicalCanonicalizationWorkspace
    search::GC.DirectedCanonicalizationWorkspace
    result::GC.DirectedCanonicalizationBuffer
    search_cache::Vector{GC.DirectedCanonicalizationWorkspace}
    result_cache::Vector{GC.DirectedCanonicalizationBuffer}
    cache_sizes::Vector{Int}
    multiplicities::Vector{Int}
    capacity::Int
end

function PhysicalCanonicalizationWorkspace(capacity::Integer)
    n = Int(capacity)
    1 <= n <= GC_PHYSICAL_MAX_VERTICES || throw(
        ArgumentError(
            "physical canonicalization workspace capacity must lie in 1:$GC_PHYSICAL_MAX_VERTICES",
        ),
    )
    search = GC.DirectedCanonicalizationWorkspace(n)
    result = GC.DirectedCanonicalizationBuffer(n)
    return PhysicalCanonicalizationWorkspace(
        search,
        result,
        GC.DirectedCanonicalizationWorkspace[search],
        GC.DirectedCanonicalizationBuffer[result],
        Int[n],
        zeros(Int, n * n),
        n,
    )
end

@inline function physical_canonicalization_workspace(::Val{E}) where {E}
    return PhysicalCanonicalizationWorkspace(
        max(1, min(GC_PHYSICAL_MAX_VERTICES, 3 * Int(E)))
    )
end

@inline function _activate_gc_workspace!(
    scratch::PhysicalCanonicalizationWorkspace, n::Int
)
    1 <= n <= scratch.capacity || throw(
        DimensionMismatch("physical canonicalization graph exceeds workspace capacity")
    )
    @inbounds for i in eachindex(scratch.cache_sizes)
        if scratch.cache_sizes[i] == n
            scratch.search = scratch.search_cache[i]
            scratch.result = scratch.result_cache[i]
            return nothing
        end
    end

    search = GC.DirectedCanonicalizationWorkspace(n)
    result = GC.DirectedCanonicalizationBuffer(n)
    push!(scratch.search_cache, search)
    push!(scratch.result_cache, result)
    push!(scratch.cache_sizes, n)
    scratch.search = search
    scratch.result = result
    return nothing
end

@inline function _clear_physical_multiplicities!(
    scratch::PhysicalCanonicalizationWorkspace, n::Int
)
    @inbounds for i in 1:(n * n)
        scratch.multiplicities[i] = 0
    end
    return nothing
end

function _gc_direct_position_graph!(
    scratch::PhysicalCanonicalizationWorkspace, vs, graph_positions::Vector{Position}
)
    n = length(graph_positions)
    _activate_gc_workspace!(scratch, n)
    _clear_physical_multiplicities!(scratch, n)
    simple = true
    @inbounds for item in vs
        out, in = positions(item)
        source = position_vertex(graph_positions, out)
        target = position_vertex(graph_positions, in)
        slot = (source - 1) * n + target
        simple &= iszero(scratch.multiplicities[slot])
        scratch.multiplicities[slot] += 1
    end
    return GC.DirectedGCGraph(n, scratch.multiplicities), position_labels(graph_positions), simple
end

function _gc_colored_subdivision_graph!(
    scratch::PhysicalCanonicalizationWorkspace, vs, graph_positions::Vector{Position}
)
    colors = propagator_colors(vs)
    npositions = length(graph_positions)
    n = npositions + length(vs)
    _activate_gc_workspace!(scratch, n)
    _clear_physical_multiplicities!(scratch, n)

    labels = Vector{Int}(undef, n)
    position_colors = position_labels(graph_positions)
    copyto!(labels, 1, position_colors, 1, npositions)

    @inbounds for (i, item) in enumerate(vs)
        color_index = searchsortedfirst(colors, propagator_color(item))
        edge_vertex = npositions + i
        labels[edge_vertex] = 3 + color_index
        out, in = positions(item)
        source = position_vertex(graph_positions, out)
        target = position_vertex(graph_positions, in)
        scratch.multiplicities[(source - 1) * n + edge_vertex] += 1
        scratch.multiplicities[(edge_vertex - 1) * n + target] += 1
    end
    return GC.DirectedGCGraph(n, scratch.multiplicities), labels
end

"""
Canonicalize the physical graph with released GraphCombinations v0.4.0.

The direct position graph is tried first. Its witness is already physically sufficient when its
exact automorphism group is trivial, or for a simple uniformly colored propagator graph. Only
residual physical ambiguity uses the colored subdivision representation. `false` means the
required representation exceeds the released one-word `n <= 64` backend and the caller must use
the existing Nauty fallback.
"""
function _gc_physical_witness!(
    scratch::PhysicalCanonicalizationWorkspace, vs, graph_positions::Vector{Position}
)::Bool
    npositions = length(graph_positions)
    npositions <= scratch.capacity || return false

    direct_graph, vertex_colors, simple =
        _gc_direct_position_graph!(scratch, vs, graph_positions)
    GC.canonicalize_directed!(scratch.result, scratch.search, direct_graph, vertex_colors)

    use_direct =
        isone(GC.canonical_automorphism_order(scratch.result)) ||
        (simple && uniform_coloring(vs))
    use_direct && return true

    ncolored = npositions + length(vs)
    ncolored <= scratch.capacity || return false
    colored_graph, colored_colors =
        _gc_colored_subdivision_graph!(scratch, vs, graph_positions)
    GC.canonicalize_directed!(scratch.result, scratch.search, colored_graph, colored_colors)
    return true
end

"""
Construct the vertex-colored directed graph used for canonicalization.

A uniform simple propagator set uses the original position graph directly: when every edge
has the same physical color and no directed position-pair is repeated, edge colors carry no
additional isomorphism information. Mixed-color and multiedge graphs use labeled subdivision
vertices so field family/index, propagator type, orientation, regularisation, and derivative
endpoint decoration remain part of the canonical form. `Out()` and `In()` always have fixed,
distinct vertex colors.
"""
function make_NautyDiGraph(vs::Vector{T}) where {T<:Union{Contraction,Edge}}
    isempty(vs) && return NautyGraphs.NautyDiGraph(0), Position[]

    graph_positions = canonicalization_positions(vs)
    simple_graph, simple = _make_simple_NautyDiGraph(vs, graph_positions)
    graph = if simple && uniform_coloring(vs)
        simple_graph
    else
        make_colored_NautyDiGraph(vs, graph_positions)
    end
    return graph, graph_positions
end
function make_NautyDiGraph(vs::Vector{Tuple{Field{S},Field{S}}}) where {S<:Statistics}
    contractions = Contraction{S}[Contraction(v) for v in vs]
    return make_NautyDiGraph(contractions)
end

"""
Return physical and uncolored canonicalization metadata from one direct-graph Nauty call.

The direct position graph defines topology. Its canonical permutation is retained independently
of physical edge colors. The same permutation is sufficient for physical canonicalization
when the direct graph has no nontrivial automorphism, or when a simple graph has uniform
physical coloring. Only genuinely ambiguous colored graphs require the subdivision fallback.
"""
function canonicalization_permutations(vs, graph_positions::Vector{Position})
    graph, simple = _make_simple_NautyDiGraph(vs, graph_positions)
    topology_permutation, automorphisms = NautyGraphs.nauty(graph)

    physical_permutation = if isone(automorphisms.n) || (simple && uniform_coloring(vs))
        topology_permutation
    else
        colored_graph = make_colored_NautyDiGraph(vs, graph_positions)
        NautyGraphs.canonical_permutation(colored_graph)
    end
    return physical_permutation, topology_permutation, simple, automorphisms
end

function canonicalization_permutation(vs, graph_positions::Vector{Position})
    physical_permutation, _, _, _ = canonicalization_permutations(vs, graph_positions)
    return physical_permutation
end

function out_bulk_positions(vs)
    result = Position[]
    for item in vs
        p1, p2 = positions(item)
        if is_out(p1) && is_bulk(p2)
            p2 in result || push!(result, p2)
        elseif is_out(p2) && is_bulk(p1)
            p1 in result || push!(result, p1)
        end
    end
    return result
end

function make_permutation_dict(
    perm::AbstractVector{<:Integer}, graph_positions::Vector{Position}, vs
)
    npositions = length(graph_positions)
    canonical_bulk = Position[]
    for original_vertex in perm
        original_vertex <= npositions || continue
        old_position = graph_positions[original_vertex]
        is_bulk(old_position) && push!(canonical_bulk, old_position)
    end

    # Preserve the package's external-anchor convention: bulk vertices attached to Out()
    # receive the first canonical labels, ordered by their backend canonical rank.
    anchors = out_bulk_positions(vs)
    mapping = Dict{Position,Position}()
    bulk_index = 0
    for old_position in canonical_bulk
        old_position in anchors || continue
        bulk_index += 1
        mapping[old_position] = Bulk(bulk_index)
    end
    for old_position in canonical_bulk
        old_position in anchors && continue
        bulk_index += 1
        mapping[old_position] = Bulk(bulk_index)
    end
    return mapping
end

function _gc_permutation_dict(
    buffer::GC.DirectedCanonicalizationBuffer, graph_positions::Vector{Position}, vs
)
    npositions = length(graph_positions)
    canonical_bulk = Position[]
    sizehint!(canonical_bulk, npositions)
    @inbounds for rank in 1:buffer.num_vertices
        original_vertex = GC.original_vertex(buffer, rank)
        original_vertex <= npositions || continue
        old_position = graph_positions[original_vertex]
        is_bulk(old_position) && push!(canonical_bulk, old_position)
    end

    anchors = out_bulk_positions(vs)
    mapping = Dict{Position,Position}()
    sizehint!(mapping, length(canonical_bulk))
    bulk_index = 0
    for old_position in canonical_bulk
        old_position in anchors || continue
        bulk_index += 1
        mapping[old_position] = Bulk(bulk_index)
    end
    for old_position in canonical_bulk
        old_position in anchors && continue
        bulk_index += 1
        mapping[old_position] = Bulk(bulk_index)
    end
    return mapping
end

function relabel_bulk_position(f::Field, mapping::Dict{Position,Position})
    p = position(f)
    return is_bulk(p) && haskey(mapping, p) ? f(mapping[p]) : f
end
function relabel_bulk_positions(
    c::Contraction{S}, mapping::Dict{Position,Position}
) where {S}
    return Contraction(
        relabel_bulk_position(c.out, mapping), relabel_bulk_position(c.in, mapping)
    )
end
function relabel_bulk_positions(e::Edge{S}, mapping::Dict{Position,Position}) where {S}
    return Edge(
        relabel_bulk_position(e.out, mapping),
        relabel_bulk_position(e.in, mapping),
        e.edgetype,
        e.momenta,
    )
end

"""
Build the uncolored direct position graph exactly as the pre-static implementation did.

This helper is intentionally separate from physical canonicalization. Its sole purpose is to
preserve the historical topology-label contract used by the analytical `1 / 3 / 11 / 59`
regressions: external vertices are encoded by their legacy graph positions, bulk vertices are
interchangeable, edge colors are ignored, and duplicate directed edges do not affect the vertex
canonical permutation.
"""
function legacy_topology_graph(vs)
    position_pairs = Tuple{Int8,Int8}[integer_positions(item) for item in vs]
    flattened = collect(Iterators.flatten(position_pairs))
    max_label = length(unique(flattened))
    has_out = typemin(Int8) in flattened

    graph_edges = map(position_pairs) do pair
        vertices = if typemin(Int8) in pair
            (1, last(pair) + Int(has_out))
        elseif typemax(Int8) in pair
            (first(pair) + Int(has_out), max_label)
        else
            pair .+ Int(has_out)
        end
        return Graphs.Edge(vertices)
    end
    return NautyGraphs.NautyDiGraph(graph_edges), max_label, has_out
end

function legacy_topology_permutation_dict(
    permutation::AbstractVector{<:Integer}, max_label::Int, has_out::Bool
)
    if has_out
        tracker = 0
        in_index = findfirst(==(max_label), permutation)
        out_index = findfirst(==(1), permutation)
        mapping = Dict{Position,Position}()
        for i in eachindex(permutation)
            if i == out_index || i == in_index
                tracker += 1
            else
                mapping[Bulk(permutation[i] - 1)] = Bulk(i - tracker)
            end
        end
        return mapping
    end

    return Dict{Position,Position}(
        Bulk(permutation[i]) => Bulk(i) for i in eachindex(permutation)
    )
end

function legacy_topology(vs, ::Val{E2}) where {E2}
    isempty(vs) && return bulk_multiplicity(Tuple{Int8,Int8}[], Val(E2))
    graph, max_label, has_out = legacy_topology_graph(vs)
    permutation = NautyGraphs.canonical_permutation(graph)
    mapping = legacy_topology_permutation_dict(permutation, max_label, has_out)
    topology_edges = Tuple{Int8,Int8}[
        integer_positions(relabel_bulk_positions(item, mapping)) for item in vs
    ]
    return bulk_multiplicity(topology_edges, Val(E2))
end

function _canonicalize_typed(
    vs::Vector{T},
    scratch::PhysicalCanonicalizationWorkspace,
    graph_positions::Vector{Position},
) where {T}
    if _gc_physical_witness!(scratch, vs, graph_positions)
        permutation_map = _gc_permutation_dict(scratch.result, graph_positions, vs)
    else
        physical_permutation = canonicalization_permutation(vs, graph_positions)
        permutation_map = make_permutation_dict(physical_permutation, graph_positions, vs)
    end
    return T[relabel_bulk_positions(item, permutation_map) for item in vs]
end

function _canonicalize_typed(
    vs::Vector{T}, scratch::PhysicalCanonicalizationWorkspace
) where {T}
    isempty(vs) && return copy(vs)
    graph_positions = canonicalization_positions(vs)
    return _canonicalize_typed(vs, scratch, graph_positions)
end

function _canonicalize_typed(vs::Vector{T}) where {T}
    isempty(vs) && return copy(vs)
    graph_positions = canonicalization_positions(vs)
    capacity = max(1, min(GC_PHYSICAL_MAX_VERTICES, length(graph_positions) + length(vs)))
    scratch = PhysicalCanonicalizationWorkspace(capacity)
    return _canonicalize_typed(vs, scratch, graph_positions)
end

canonicalize(vs::Vector{Contraction{S}}) where {S<:Statistics} = _canonicalize_typed(vs)
canonicalize(vs::Vector{Edge{S}}) where {S<:Statistics} = _canonicalize_typed(vs)
canonicalize(
    vs::Vector{Contraction{S}}, scratch::PhysicalCanonicalizationWorkspace
) where {S<:Statistics} = _canonicalize_typed(vs, scratch)
canonicalize(
    vs::Vector{Edge{S}}, scratch::PhysicalCanonicalizationWorkspace
) where {S<:Statistics} = _canonicalize_typed(vs, scratch)

"""
Canonicalize physical contractions and compute the established uncolored topology signature.

Physical canonicalization uses released GraphCombinations v0.4.0 for every supported direct or
colored-subdivision graph and falls back to the existing Nauty representation only beyond the
one-word `n <= 64` backend. Topology metadata remains independently computed with the historical
Nauty path, so physical field colors cannot fragment topology classes and the established topology
labels remain exact.
"""
function _canonicalize_with_topology_typed(vs::Vector{T}, ::Val{E2}) where {T,E2}
    isempty(vs) && return copy(vs), bulk_multiplicity(Tuple{Int8,Int8}[], Val(E2))

    canonical_vs = _canonicalize_typed(vs)
    topology = legacy_topology(vs, Val(E2))
    return canonical_vs, topology
end

function canonicalize_with_topology(
    vs::Vector{Contraction{S}}, ::Val{E2}
) where {S<:Statistics,E2}
    return _canonicalize_with_topology_typed(vs, Val(E2))
end
function canonicalize_with_topology(vs::Vector{Edge{S}}, ::Val{E2}) where {S<:Statistics,E2}
    return _canonicalize_with_topology_typed(vs, Val(E2))
end

function canonicalize(vs::Vector{Tuple{Field{S},Field{S}}}) where {S<:Statistics}
    contractions = Contraction{S}[Contraction(v) for v in vs]
    return Tuple{Field{S},Field{S}}[Tuple(c) for c in canonicalize(contractions)]
end