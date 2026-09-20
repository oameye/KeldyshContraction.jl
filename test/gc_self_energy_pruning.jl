using KeldyshContraction, Test
import KeldyshContraction as KC
const GraphComb1PI = KC.GC

function _oracle_connected(edges::Vector{Tuple{Int,Int}})::Bool
    isempty(edges) && return true
    vertices = Set{Int}()
    for edge in edges
        push!(vertices, edge[1])
        push!(vertices, edge[2])
    end
    length(vertices) <= 1 && return true

    source = first(vertices)
    visited = Set{Int}((source,))
    queue = Int[source]
    while !isempty(queue)
        current = popfirst!(queue)
        for edge in edges
            neighbor = if edge[1] == current
                edge[2]
            elseif edge[2] == current
                edge[1]
            else
                0
            end
            iszero(neighbor) && continue
            if neighbor ∉ visited
                push!(visited, neighbor)
                push!(queue, neighbor)
            end
        end
    end
    return length(visited) == length(vertices)
end

function _oracle_irreducible(edges::Vector{Tuple{Int,Int}})::Bool
    length(edges) < 2 && return true
    for skipped in eachindex(edges)
        remaining = Tuple{Int,Int}[edges[i] for i in eachindex(edges) if i != skipped]
        _oracle_connected(remaining) || return false
    end
    return true
end

function _has_irreducible_completion(
    policy::KC._GCOnePIPruningPolicy,
    edges::Vector{GraphComb1PI.ColoredPortEdge},
    sources::Matrix{Int},
    targets::Matrix{Int},
)::Bool
    source_cell = nothing
    for vertex in axes(sources, 1), color in axes(sources, 2)
        if !iszero(sources[vertex, color])
            source_cell = (vertex, color)
            break
        end
    end

    if isnothing(source_cell)
        all(iszero, targets) || return false
        bulk_edges = Tuple{Int,Int}[]
        for edge in edges
            policy.bulk_allowed[
                edge.source, edge.source_color, edge.target, edge.target_color
            ] || continue
            push!(bulk_edges, (edge.source, edge.target))
        end
        return _oracle_irreducible(bulk_edges)
    end

    source_vertex, source_color = source_cell
    next_sources = copy(sources)
    next_sources[source_vertex, source_color] -= 1
    for target_vertex in axes(targets, 1), target_color in axes(targets, 2)
        iszero(targets[target_vertex, target_color]) && continue
        next_targets = copy(targets)
        next_targets[target_vertex, target_color] -= 1
        edge = GraphComb1PI.ColoredPortEdge(
            source_vertex, target_vertex, source_color, target_color
        )
        push!(edges, edge)
        if _has_irreducible_completion(policy, edges, next_sources, next_targets)
            pop!(edges)
            return true
        end
        pop!(edges)
    end
    return false
end

function _certify_reachable_states!(
    policy::KC._GCOnePIPruningPolicy,
    edges::Vector{GraphComb1PI.ColoredPortEdge},
    sources::Matrix{Int},
    targets::Matrix{Int},
    checked::Base.RefValue{Int},
)
    state = GraphComb1PI.ColoredPortState(edges, sources, targets)
    checked[] += 1
    if !policy(state)
        @test !_has_irreducible_completion(
            policy, copy(edges), copy(sources), copy(targets)
        )
        return nothing
    end

    source_cell = nothing
    for vertex in axes(sources, 1), color in axes(sources, 2)
        if !iszero(sources[vertex, color])
            source_cell = (vertex, color)
            break
        end
    end
    isnothing(source_cell) && return nothing

    source_vertex, source_color = source_cell
    next_sources = copy(sources)
    next_sources[source_vertex, source_color] -= 1
    for target_vertex in axes(targets, 1), target_color in axes(targets, 2)
        iszero(targets[target_vertex, target_color]) && continue
        next_targets = copy(targets)
        next_targets[target_vertex, target_color] -= 1
        edge = GraphComb1PI.ColoredPortEdge(
            source_vertex, target_vertex, source_color, target_color
        )
        push!(edges, edge)
        _certify_reachable_states!(policy, edges, next_sources, next_targets, checked)
        pop!(edges)
    end
end

function _signed_pairing_weights(pairings)
    weights = Dict{Any,Int}()
    for (pairing, topology, multiplicity) in pairings
        key = (pairing.contractions, topology)
        weight = Int(pairing.sign) * multiplicity
        weights[key] = get(weights, key, 0) + weight
    end
    filter!(pair -> !iszero(last(pair)), weights)
    return weights
end

@testset "Continuation-safe 1PI pruning" begin
    # KC's oracle drops isolated vertices after edge removal. A two-edge chain is therefore
    # irreducible by the existing package semantics and must not be pruned.
    @test !KC._has_unhealable_bulk_bridge([(1, 2), (2, 3)], Tuple{Int,Int}[])

    # The middle edge of a three-edge chain separates two persistent edge-containing components.
    @test KC._has_unhealable_bulk_bridge([(1, 2), (2, 3), (3, 4)], Tuple{Int,Int}[])

    # A possible future edge can heal that cut, so current reducibility is not monotone.
    @test !KC._has_unhealable_bulk_bridge([(1, 2), (2, 3), (3, 4)], [(1, 4)])

    @test !KC._has_unhealable_bulk_bridge([(1, 2), (2, 3), (3, 1)], Tuple{Int,Int}[])
    @test !KC._has_unhealable_bulk_bridge([(1, 2), (1, 2), (2, 3)], Tuple{Int,Int}[])

    # Exhaust all 3x3 one-color bulk/nonbulk classifications. Every port pairing is allowed;
    # a false bulk cell therefore models a nonbulk contraction that may consume residual ports but
    # cannot heal a 1PI cut. This is stricter than treating false cells as simply forbidden.
    # Keep this oracle eager at every layer; production depth gating is tested separately below.
    checked = Ref(0)
    source_ports = ones(Int, 3, 1)
    target_ports = ones(Int, 3, 1)
    for mask in 0:(2 ^ 9 - 1)
        bulk_allowed = falses(3, 1, 3, 1)
        bit = 0
        for source in 1:3, target in 1:3
            bulk_allowed[source, 1, target, 1] = !iszero(mask & (1 << bit))
            bit += 1
        end
        policy = KC._GCOnePIPruningPolicy(bulk_allowed, 3, 3)
        _certify_reachable_states!(
            policy, GraphComb1PI.ColoredPortEdge[], source_ports, target_ports, checked
        )
    end
    @test checked[] > 5_000

    # Repeated residual ports exercise parallel propagators and multiplicity-shaped states.
    repeated_policy = KC._GCOnePIPruningPolicy(trues(2, 1, 2, 1), 3, 3)
    repeated_checked = Ref(0)
    _certify_reachable_states!(
        repeated_policy,
        GraphComb1PI.ColoredPortEdge[],
        reshape([2, 1], 2, 1),
        reshape([2, 1], 2, 1),
        repeated_checked,
    )
    @test repeated_checked[] > 1

    # The opt-in physical path must equal ordinary GC Wick generation followed by the established
    # irreducibility oracle. This does not alter production routing; `onepi_pruning` defaults false.
    @qfields gc1pi_ϕ::Boson
    c, q = gc1pi_ϕ[Classical], gc1pi_ϕ[Quantum]
    elastic = -(0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2))
    L = InteractionLagrangian(elastic)
    expression = c(KC.Out()) * bar(c)(KC.In()) * L(1).lagrangian * L(2).lagrangian

    for term in KC.terms(expression)
        ordinary, _ = KC._gc_wick_contraction_with_stats(
            term.args_nc, Val(5), Val(1); regularise=true, simplify=false
        )
        expected = [
            entry for entry in ordinary if KC.is_irreducible(first(entry).contractions)
        ]
        direct_onepi, _ = KC._gc_wick_contraction_with_stats(
            term.args_nc,
            Val(5),
            Val(1);
            regularise=true,
            simplify=false,
            onepi_pruning=true,
        )
        @test _signed_pairing_weights(direct_onepi) == _signed_pairing_weights(expected)
    end
end
