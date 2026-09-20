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
"""
function _has_unhealable_bulk_bridge(
    existing_edges::Vector{Tuple{Int,Int}}, possible_edges::Vector{Tuple{Int,Int}}
)::Bool
    @inbounds for i in eachindex(existing_edges)
        source, target = existing_edges[i]
        source == target && continue
        _bulk_path_exists_without_edge(
            existing_edges, possible_edges, i, source, target
        ) || return true
    end
    return false
end
