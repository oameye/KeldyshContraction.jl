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

@inline function field_color(f::Field)
    return (
        name(f),
        slots(field_indices(f)),
        Int(orientation(f)),
        Int(keldysh_index(f)),
        Int(regularisation(f)),
    )
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

function uniform_coloring(vs)
    isempty(vs) && return true
    color = propagator_color(first(vs))
    return all(item -> isequal(propagator_color(item), color), vs)
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

"""
Construct the vertex-colored directed graph used for canonicalization.

A uniform simple propagator set uses the original position graph directly: when every edge
has the same physical color and no directed position-pair is repeated, edge colors carry no
additional isomorphism information. Mixed-color and multiedge graphs use labeled subdivision
vertices so field family/index, propagator type, orientation, and regularisation remain part
of the canonical form. `Out()` and `In()` always have fixed, distinct vertex colors.
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
    # receive the first canonical labels, ordered by their Nauty canonical rank.
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

function canonicalize(vs::Vector{T}) where {T<:Union{Contraction,Edge}}
    isempty(vs) && return copy(vs)
    graph_positions = canonicalization_positions(vs)
    physical_permutation, _, _, _ = canonicalization_permutations(vs, graph_positions)
    permutation_map = make_permutation_dict(physical_permutation, graph_positions, vs)
    return T[relabel_bulk_positions(item, permutation_map) for item in vs]
end

"""
Canonicalize physical contractions and compute the established uncolored topology signature.

Topology is defined after physical canonical labels have been fixed. For simple position graphs,
or for multigraphs with a trivial uncolored automorphism group, the resulting multiplicity
signature can be read directly from the physically canonicalized contractions. Only a symmetric
multigraph needs a second direct-graph canonicalization to reproduce the historical topology
tie breaking. This preserves the analytical `1 / 3 / 11 / 59` topology contract without paying
a redundant Nauty pass for every Wick matching.
"""
function canonicalize_with_topology(
    vs::Vector{T}, ::Val{E2}
) where {T<:Union{Contraction,Edge},E2}
    isempty(vs) && return copy(vs), bulk_multiplicity(Tuple{Int8,Int8}[], Val(E2))

    graph_positions = canonicalization_positions(vs)
    physical_permutation, _, simple, automorphisms =
        canonicalization_permutations(vs, graph_positions)
    physical_map = make_permutation_dict(physical_permutation, graph_positions, vs)
    canonical_vs = T[relabel_bulk_positions(item, physical_map) for item in vs]

    topology_edges = if simple || isone(automorphisms.n)
        Tuple{Int8,Int8}[integer_positions(item) for item in canonical_vs]
    else
        canonical_positions = canonicalization_positions(canonical_vs)
        graph = make_simple_NautyDiGraph(canonical_vs, canonical_positions)
        topology_permutation = NautyGraphs.canonical_permutation(graph)
        topology_map =
            make_permutation_dict(topology_permutation, canonical_positions, canonical_vs)
        Tuple{Int8,Int8}[
            integer_positions(relabel_bulk_positions(item, topology_map)) for item in canonical_vs
        ]
    end

    topology = bulk_multiplicity(topology_edges, Val(E2))
    return canonical_vs, topology
end

function canonicalize(vs::Vector{Tuple{Field{S},Field{S}}}) where {S<:Statistics}
    contractions = Contraction{S}[Contraction(v) for v in vs]
    return Tuple{Field{S},Field{S}}[Tuple(c) for c in canonicalize(contractions)]
end
