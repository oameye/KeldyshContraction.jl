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
    total_edges::Int
    max_residual_pairs::Int
end

struct _GCAndPortPolicy{A,B}
    first::A
    second::B
end

@inline function (policy::_GCAndPortPolicy)(state::GC.ColoredPortState)::Bool
    return policy.first(state) && policy.second(state)
end

function _gc_onepi_policy(
    lookup::Dict{NTuple{4,Int},Contraction{S}}, ::Val{E}; max_residual_pairs=1
) where {S<:Statistics,E}
    isempty(lookup) &&
        return _GCOnePIPruningPolicy(falses(0, 0, 0, 0), E, max_residual_pairs)

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
    return _GCOnePIPruningPolicy(bulk_allowed, E, max_residual_pairs)
end

@inline function _gc_bulk_allowed(
    policy::_GCOnePIPruningPolicy,
    source_vertex::Int,
    source_color::Int,
    target_vertex::Int,
    target_color::Int,
)::Bool
    allowed = policy.bulk_allowed
    source_vertex in axes(allowed, 1) || return false
    source_color in axes(allowed, 2) || return false
    target_vertex in axes(allowed, 3) || return false
    target_color in axes(allowed, 4) || return false
    return allowed[source_vertex, source_color, target_vertex, target_color]
end

@inline function _gc_bulk_cell(
    policy::_GCOnePIPruningPolicy, edge::GC.ColoredPortEdge
)::Bool
    return _gc_bulk_allowed(
        policy, edge.source, edge.source_color, edge.target, edge.target_color
    )
end

@inline function _gc_connect_onepi_vertices(
    reachable::UInt128, source::Int, target::Int
)::UInt128
    source_bit = _gc_port_vertex_bit(source)
    target_bit = _gc_port_vertex_bit(target)
    if !iszero(reachable & (source_bit | target_bit))
        return reachable | source_bit | target_bit
    end
    return reachable
end

@inline function _gc_single_residual_bulk_edge(
    policy::_GCOnePIPruningPolicy, state::GC.ColoredPortState
)::Tuple{Int,Int}
    sources = GC.source_port_counts(state)
    targets = GC.target_port_counts(state)
    source_vertex = 0
    source_color = 0
    target_vertex = 0
    target_color = 0

    for vertex in axes(sources, 1), color in axes(sources, 2)
        iszero(sources[vertex, color]) && continue
        source_vertex = vertex
        source_color = color
        break
    end
    iszero(source_vertex) && return 0, 0

    for vertex in axes(targets, 1), color in axes(targets, 2)
        iszero(targets[vertex, color]) && continue
        target_vertex = vertex
        target_color = color
        break
    end
    iszero(target_vertex) && return 0, 0

    _gc_bulk_allowed(policy, source_vertex, source_color, target_vertex, target_color) ||
        return 0, 0
    return source_vertex, target_vertex
end

function _gc_onepi_reachable_single_residual(
    policy::_GCOnePIPruningPolicy,
    state::GC.ColoredPortState,
    skipped_bulk_edge::Int,
    source::Int,
    possible_source::Int,
    possible_target::Int,
)::UInt128
    edges = GC.port_edges(state)
    reachable = _gc_port_vertex_bit(source)

    previous = zero(UInt128)
    while reachable != previous
        previous = reachable

        bulk_index = 0
        for edge in edges
            _gc_bulk_cell(policy, edge) || continue
            bulk_index += 1
            bulk_index == skipped_bulk_edge && continue
            reachable = _gc_connect_onepi_vertices(reachable, edge.source, edge.target)
        end

        iszero(possible_source) ||
            (reachable = _gc_connect_onepi_vertices(reachable, possible_source, possible_target))
    end
    return reachable
end

function _gc_onepi_reachable_generic(
    policy::_GCOnePIPruningPolicy,
    state::GC.ColoredPortState,
    skipped_bulk_edge::Int,
    source::Int,
)::UInt128
    edges = GC.port_edges(state)
    sources = GC.source_port_counts(state)
    targets = GC.target_port_counts(state)
    reachable = _gc_port_vertex_bit(source)

    previous = zero(UInt128)
    while reachable != previous
        previous = reachable

        bulk_index = 0
        for edge in edges
            _gc_bulk_cell(policy, edge) || continue
            bulk_index += 1
            bulk_index == skipped_bulk_edge && continue
            reachable = _gc_connect_onepi_vertices(reachable, edge.source, edge.target)
        end

        for source_vertex in axes(sources, 1)
            for source_color in axes(sources, 2)
                iszero(sources[source_vertex, source_color]) && continue
                for target_vertex in axes(targets, 1)
                    for target_color in axes(targets, 2)
                        iszero(targets[target_vertex, target_color]) && continue
                        _gc_bulk_allowed(
                            policy, source_vertex, source_color, target_vertex, target_color
                        ) || continue
                        reachable = _gc_connect_onepi_vertices(
                            reachable, source_vertex, target_vertex
                        )
                    end
                end
            end
        end
    end
    return reachable
end

function (policy::_GCOnePIPruningPolicy)(state::GC.ColoredPortState)::Bool
    edges = GC.port_edges(state)
    residual_pairs = policy.total_edges - length(edges)
    residual_pairs > policy.max_residual_pairs && return true

    bulk_count = 0
    for edge in edges
        _gc_bulk_cell(policy, edge) && (bulk_count += 1)
    end
    bulk_count < 3 && return true

    possible_source, possible_target = if residual_pairs <= 1
        _gc_single_residual_bulk_edge(policy, state)
    else
        (0, 0)
    end

    for skipped in 1:bulk_count
        bulk_index = 0
        source = 0
        first_remaining = 0
        for edge in edges
            _gc_bulk_cell(policy, edge) || continue
            bulk_index += 1
            bulk_index == skipped && continue
            source = edge.source
            first_remaining = bulk_index
            break
        end

        iszero(source) && continue
        reachable = if residual_pairs <= 1
            _gc_onepi_reachable_single_residual(
                policy, state, skipped, source, possible_source, possible_target
            )
        else
            _gc_onepi_reachable_generic(policy, state, skipped, source)
        end
        bulk_index = 0
        for edge in edges
            _gc_bulk_cell(policy, edge) || continue
            bulk_index += 1
            (bulk_index == skipped || bulk_index == first_remaining) && continue
            iszero(reachable & _gc_port_vertex_bit(edge.source)) && return false
        end
    end
    return true
end
