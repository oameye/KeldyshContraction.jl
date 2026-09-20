# KC-owned continuation-safe structural pruning for direct self-energy generation.
#
# KC's current `is_irreducible` semantics remove one bulk contraction and then recompute graph
# vertices from the remaining edges. Consequently a graph-theoretic bridge is not by itself a safe
# pruning witness: if removing it leaves only one edge-containing component, any isolated endpoint
# disappears from the oracle graph. We may reject only when a realized edge has an unavoidable cut
# separating at least two persistent edge-containing components under every continuation.

@inline function _other_bridge_vertex(edge::Tuple{Int,Int}, vertex::Int)
    if edge[1] == vertex
        return edge[2]
    elseif edge[2] == vertex
        return edge[1]
    end
    return 0
end

function _bulk_reachable_vertices(
    existing_edges::Vector{Tuple{Int,Int}},
    possible_edges::Vector{Tuple{Int,Int}},
    skipped_edge::Int,
    source::Int,
)::Set{Int}
    visited = Set{Int}((source,))
    queue = Int[source]
    while !isempty(queue)
        current = popfirst!(queue)

        @inbounds for i in eachindex(existing_edges)
            i == skipped_edge && continue
            neighbor = _other_bridge_vertex(existing_edges[i], current)
            iszero(neighbor) && continue
            if neighbor ∉ visited
                push!(visited, neighbor)
                push!(queue, neighbor)
            end
        end

        for edge in possible_edges
            neighbor = _other_bridge_vertex(edge, current)
            iszero(neighbor) && continue
            if neighbor ∉ visited
                push!(visited, neighbor)
                push!(queue, neighbor)
            end
        end
    end
    return visited
end

function _has_unhealable_bulk_bridge(
    existing_edges::Vector{Tuple{Int,Int}}, possible_edges::Vector{Tuple{Int,Int}}
)::Bool
    # After removing one edge, KC's oracle can only observe disconnectedness if at least two
    # already-realized bulk edges remain. Future edges are optional in the optimistic support and
    # therefore cannot be used to establish persistence of a component.
    length(existing_edges) < 3 && return false

    @inbounds for skipped in eachindex(existing_edges)
        first_remaining = findfirst(i -> i != skipped, eachindex(existing_edges))
        isnothing(first_remaining) && continue
        source = existing_edges[first_remaining][1]
        reachable = _bulk_reachable_vertices(
            existing_edges, possible_edges, skipped, source
        )

        for i in eachindex(existing_edges)
            (i == skipped || i == first_remaining) && continue
            edge = existing_edges[i]
            if edge[1] ∉ reachable
                return true
            end
        end
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
    length(existing) < 3 && return true
    possible = _gc_possible_bulk_edges(policy, state)
    return !_has_unhealable_bulk_bridge(existing, possible)
end
