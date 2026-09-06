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
function sort_by_position_and_type(p::Tuple{Field{S},Field{S}})::Float64 where {S<:Statistics}
    return sort_by_position_and_type(Contraction(p))
end
sort_by_position_and_type(p::Edge)::Float64 =
    sort_by_position_and_type(Contraction(fields(p)))

field_color(f::Field) = (
    string(name(f)),
    slots(field_indices(f)),
    Int(orientation(f)),
    Int(keldysh_index(f)),
    Int(regularisation(f)),
)

function propagator_color(c::Contraction)
    return (field_color(c.out), field_color(c.in), Int(propagator_type(c...)))
end
function propagator_color(e::Edge)
    return (field_color(e.out), field_color(e.in), Int(propagator_type(e)))
end

"""
Construct the vertex-colored directed graph used for canonicalization.

Physical propagator colors are represented by labeled subdivision vertices because
NautyGraphs does not support edge labels. `Out()` and `In()` receive distinct vertex
colors, while all bulk integration vertices share one color and may be relabeled.
"""
function make_NautyDiGraph(vs::Vector{T}) where {T<:Union{Contraction,Edge}}
    isempty(vs) && return NautyGraphs.NautyDiGraph(0), Position[]

    graph_positions = sort!(unique(Position[p for item in vs for p in positions(item)]))
    position_vertex = Dict(p => i for (i, p) in enumerate(graph_positions))

    colors = sort!(unique(propagator_color.(vs)))
    npositions = length(graph_positions)
    labels = Vector{Int}(undef, npositions + length(vs))

    for (i, p) in enumerate(graph_positions)
        labels[i] = is_out(p) ? 1 : is_in(p) ? 2 : 3
    end
    for (i, item) in enumerate(vs)
        color = propagator_color(item)
        color_index = findfirst(isequal(color), colors)::Int
        labels[npositions + i] = 3 + color_index
    end

    graph = NautyGraphs.NautyDiGraph(length(labels); vertex_labels=labels)
    for (i, item) in enumerate(vs)
        out, in = positions(item)
        edge_vertex = npositions + i
        Graphs.add_edge!(graph, position_vertex[out], edge_vertex)
        Graphs.add_edge!(graph, edge_vertex, position_vertex[in])
    end
    return graph, graph_positions
end
function make_NautyDiGraph(vs::Vector{Tuple{Field{S},Field{S}}}) where {S<:Statistics}
    contractions = Contraction{S}[Contraction(v) for v in vs]
    return make_NautyDiGraph(contractions)
end

function make_permutation_dict(perm::AbstractVector{<:Integer}, graph_positions::Vector{Position})
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
function relabel_bulk_positions(c::Contraction{S}, mapping::Dict{Position,Position}) where {S}
    return Contraction(relabel_bulk_position(c.out, mapping), relabel_bulk_position(c.in, mapping))
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
