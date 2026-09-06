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

function uniform_simple_coloring(vs)
    isempty(vs) && return true
    color = propagator_color(first(vs))
    for i in eachindex(vs)
        item = vs[i]
        isequal(propagator_color(item), color) || return false
        item_positions = positions(item)
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

function make_simple_NautyDiGraph(vs, graph_positions::Vector{Position})
    graph = NautyGraphs.NautyDiGraph(
        length(graph_positions); vertex_labels=position_labels(graph_positions)
    )
    for item in vs
        out, in = positions(item)
        Graphs.add_edge!(
            graph,
            position_vertex(graph_positions, out),
            position_vertex(graph_positions, in),
        )
    end
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
    graph = if uniform_simple_coloring(vs)
        make_simple_NautyDiGraph(vs, graph_positions)
    else
        make_colored_NautyDiGraph(vs, graph_positions)
    end
    return graph, graph_positions
end
function make_NautyDiGraph(vs::Vector{Tuple{Field{S},Field{S}}}) where {S<:Statistics}
    contractions = Contraction{S}[Contraction(v) for v in vs]
    return make_NautyDiGraph(contractions)
end

function make_permutation_dict(
    perm::AbstractVector{<:Integer}, graph_positions::Vector{Position}
)
    npositions = length(graph_positions)
    mapping = Dict{Position,Position}()
    bulk_index = 0
    for original_vertex in perm
        original_vertex <= npositions || continue
        old_position = graph_positions[original_vertex]
        is_bulk(old_position) || continue
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
    graph, graph_positions = make_NautyDiGraph(vs)
    perm = NautyGraphs.canonical_permutation(graph)
    permutation_map = make_permutation_dict(perm, graph_positions)
    return T[relabel_bulk_positions(item, permutation_map) for item in vs]
end
function canonicalize(vs::Vector{Tuple{Field{S},Field{S}}}) where {S<:Statistics}
    contractions = Contraction{S}[Contraction(v) for v in vs]
    return Tuple{Field{S},Field{S}}[Tuple(c) for c in canonicalize(contractions)]
end
