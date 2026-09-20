# KC-owned continuation-safe structural pruning for direct self-energy generation.
#
# A currently realized bulk bridge is not sufficient for pruning: a later Wick contraction may
# connect its two sides. It is safe to reject only when the bridge remains unavoidable even after
# adding every still-admissible residual bulk connection.

@inline function _other_bridge_vertex(edge::Tuple{Int,Int}, vertex::Int)
    if edge[1] == vertex
        return edge[2]
    elseif edge[2] == vertex
        return edge[1]
    end
    return 0
end

function _bulk_path_exists_without_edge(
    existing_edges::Vector{Tuple{Int,Int}},
    possible_edges::Vector{Tuple{Int,Int}},
    skipped_edge::Int,
    source::Int,
    target::Int,
)::Bool
    source == target && return true

    visited = Set{Int}((source,))
    queue = Int[source]
    while !isempty(queue)
        current = popfirst!(queue)

        @inbounds for i in eachindex(existing_edges)
            i == skipped_edge && continue
            neighbor = _other_bridge_vertex(existing_edges[i], current)
            iszero(neighbor) && continue
            neighbor == target && return true
            if neighbor ∉ visited
                push!(visited, neighbor)
                push!(queue, neighbor)
            end
        end

        for edge in possible_edges
            neighbor = _other_bridge_vertex(edge, current)
            iszero(neighbor) && continue
            neighbor == target && return true
            if neighbor ∉ visited
                push!(visited, neighbor)
                push!(queue, neighbor)
            end
        end
    end
    return false
end

"""
    _has_unhealable_bulk_bridge(existing_edges, possible_edges)

Return whether a bulk edge already present in a partial Wick state must remain a bridge under
all admissible continuations represented by `possible_edges`.

`possible_edges` is an optimistic support graph: it may contain mutually incompatible future
contractions. This can only make the test weaker. Therefore a `true` result is continuation-safe:
if an existing edge is still a bridge after every possible residual bulk connection is added,
no physical completion can become one-particle irreducible.

The existing `SelfEnergy` oracle treats fewer than two bulk contractions as irreducible, so the
partial-state predicate deliberately does not reject until at least two bulk edges are realized.
"""
function _has_unhealable_bulk_bridge(
    existing_edges::Vector{Tuple{Int,Int}}, possible_edges::Vector{Tuple{Int,Int}}
)::Bool
    length(existing_edges) < 2 && return false

    @inbounds for i in eachindex(existing_edges)
        source, target = existing_edges[i]
        source == target && continue
        _bulk_path_exists_without_edge(
            existing_edges, possible_edges, i, source, target
        ) || return true
    end
    return false
end

struct _GCOnePIPruningPolicy
    bulk_allowed::Array{Bool,4}
end

function _gc_onepi_policy(
    lookup::Dict{NTuple{4,Int},Contraction{S}}
) where {S<:Statistics}
    isempty(lookup) && return _GCOnePIPruningPolicy(falses(0, 0, 0, 0))

    max_source_vertex = maximum(cell[1] for cell in keys(lookup))
    max_source_color = maximum(cell[2] for cell in keys(lookup))
    max_target_vertex = maximum(cell[3] for cell in keys(lookup))
    max_target_color = maximum(cell[4] for cell in keys(lookup))
    bulk_allowed = falses(
        max_source_vertex, max_source_color, max_target_vertex, max_target_color
    )
    for (cell, contraction) in lookup
        bulk_allowed[cell...] = is_bulk(contraction)
    end
    return _GCOnePIPruningPolicy(bulk_allowed)
end

@inline function _gc_bulk_cell(
    policy::_GCOnePIPruningPolicy, edge::GC.ColoredPortEdge
)::Bool
    return policy.bulk_allowed[
        edge.source, edge.source_color, edge.target, edge.target_color
    ]
end

function _gc_existing_bulk_edges(
    policy::_GCOnePIPruningPolicy, state::GC.ColoredPortState
)::Vector{Tuple{Int,Int}}
    result = Tuple{Int,Int}[]
    for edge in GC.port_edges(state)
        _gc_bulk_cell(policy, edge) || continue
        push!(result, (edge.source, edge.target))
    end
    return result
end

function _gc_possible_bulk_edges(
    policy::_GCOnePIPruningPolicy, state::GC.ColoredPortState
)::Vector{Tuple{Int,Int}}
    sources = GC.source_port_counts(state)
    targets = GC.target_port_counts(state)
    possible = Tuple{Int,Int}[]

    for source_vertex in axes(sources, 1)
        for source_color in axes(sources, 2)
            iszero(sources[source_vertex, source_color]) && continue
            for target_vertex in axes(targets, 1)
                for target_color in axes(targets, 2)
                    iszero(targets[target_vertex, target_color]) && continue
                    policy.bulk_allowed[
                        source_vertex, source_color, target_vertex, target_color
                    ] || continue
                    edge = (source_vertex, target_vertex)
                    edge in possible || push!(possible, edge)
                end
            end
        end
    end
    return possible
end

function (policy::_GCOnePIPruningPolicy)(state::GC.ColoredPortState)::Bool
    existing = _gc_existing_bulk_edges(policy, state)
    length(existing) < 2 && return true
    possible = _gc_possible_bulk_edges(policy, state)
    return !_has_unhealable_bulk_bridge(existing, possible)
end
