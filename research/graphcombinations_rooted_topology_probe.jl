using KeldyshContraction
using Test

import GraphCombinations as GC
import KeldyshContraction as KC

const GC_SHA = "e0e5d7dc4e14f4fc356ffeee3809fdb5b47fcf12"

function gc_result(edges::Vector{Tuple{Int,Int}}, num_vertices::Int)
    graph = GC.DirectedGCGraph(
        Pair{Int,Int}[source => target for (source, target) in edges], num_vertices
    )
    return graph, GC.canonicalize_directed(graph)
end

function legacy_graph_data(vs; preserve_multiplicity::Bool)
    position_pairs = Tuple{Int8,Int8}[KC.integer_positions(item) for item in vs]
    flattened = collect(Iterators.flatten(position_pairs))
    max_label = length(unique(flattened))
    has_out = typemin(Int8) in flattened

    edges = Tuple{Int,Int}[]
    sizehint!(edges, length(position_pairs))
    for pair in position_pairs
        vertices = if typemin(Int8) in pair
            (1, Int(last(pair)) + Int(has_out))
        elseif typemax(Int8) in pair
            (Int(first(pair)) + Int(has_out), max_label)
        else
            (Int(first(pair)) + Int(has_out), Int(last(pair)) + Int(has_out))
        end
        if preserve_multiplicity || !(vertices in edges)
            push!(edges, vertices)
        end
    end
    return edges, max_label, has_out
end

function gc_legacy_topology(vs, ::Val{E2}; preserve_multiplicity::Bool) where {E2}
    isempty(vs) && return (
        topology=Tuple(KC.bulk_multiplicity(Tuple{Int8,Int8}[], Val(E2))),
        canonical=(),
        edges=Tuple{Int,Int}[],
    )

    edges, max_label, has_out = legacy_graph_data(vs; preserve_multiplicity)
    graph, result = gc_result(edges, max_label)
    old_to_canonical = GC.vertex_mapping(GC.canonical_relabeling(result))
    canonical_to_old = sortperm(old_to_canonical)
    mapping = KC.legacy_topology_permutation_dict(canonical_to_old, max_label, has_out)
    topology_edges = Tuple{Int8,Int8}[
        KC.integer_positions(KC.relabel_bulk_positions(item, mapping)) for item in vs
    ]
    topology = Tuple(KC.bulk_multiplicity(topology_edges, Val(E2)))
    canonical = Tuple(GC.canonical_graph(result).multiplicities)

    buffer = GC.DirectedCanonicalizationBuffer(graph.num_vertices)
    workspace = GC.DirectedCanonicalizationWorkspace(graph.num_vertices)
    GC.canonicalize_directed!(buffer, workspace, graph)
    @test Tuple(buffer.canonical_multiplicities) == canonical
    @test buffer.old_to_canonical == old_to_canonical

    return (topology=topology, canonical=canonical, edges=edges)
end

@qfields rooted_probe_ϕ::Boson
c, q = rooted_probe_ϕ[Classical], rooted_probe_ϕ[Quantum]

println("GC rooted topology probe ($GC_SHA)")
flush(stdout)

@testset "rooted GC topology partition versus historical Nauty partition" begin
    elastic = -(
        1 // 2 * (c^2 + q^2) * bar(c) * bar(q) + 1 // 2 * c * q * (bar(c)^2 + bar(q)^2)
    )
    L = InteractionLagrangian(elastic)
    cases = ((1, 3, 1), (2, 5, 3), (3, 7, 11), (4, 9, 59))

    for (order, edge_count, expected) in cases
        G = DressedPropagator(L, Val(order), Val(edge_count))
        component = KC.topologies(G.keldysh)
        @test length(keys(component)) == expected

        simple_to_nauty = Dict{Any,Any}()
        multigraph_to_nauty = Dict{Any,Any}()
        nauty_to_multigraph = Dict{Any,Any}()
        first_simple_split = nothing
        first_simple_collapse = nothing

        for (key, diagrams) in component
            nauty_key = Tuple(key)
            simple_keys = Set{Any}()
            multigraph_keys = Set{Any}()
            first_simple = nothing

            for diagram in diagrams
                contractions = KC.Contraction{Boson}[
                    (edge.out, edge.in) for edge in KC.contractions(diagram)
                ]
                simple = gc_legacy_topology(
                    contractions, Val(length(key)); preserve_multiplicity=false
                )
                multigraph = gc_legacy_topology(
                    contractions, Val(length(key)); preserve_multiplicity=true
                )
                push!(simple_keys, simple.topology)
                push!(multigraph_keys, multigraph.topology)
                first_simple === nothing && (first_simple = simple)
            end

            if length(simple_keys) != 1 && first_simple_split === nothing
                first_simple_split = (nauty_key=nauty_key, simple_keys=collect(simple_keys))
            end
            @test length(multigraph_keys) == 1
            multigraph_key = only(multigraph_keys)
            @test get!(nauty_to_multigraph, nauty_key, multigraph_key) == multigraph_key
            previous = get(multigraph_to_nauty, multigraph_key, nothing)
            @test previous === nothing || previous == nauty_key
            multigraph_to_nauty[multigraph_key] = nauty_key

            if length(simple_keys) == 1
                simple_key = only(simple_keys)
                previous_simple = get(simple_to_nauty, simple_key, nothing)
                if previous_simple !== nothing &&
                    previous_simple != nauty_key &&
                    first_simple_collapse === nothing
                    first_simple_collapse = (
                        first_nauty_key=previous_simple,
                        second_nauty_key=nauty_key,
                        simple_key=simple_key,
                        canonical=first_simple.canonical,
                        edges=first_simple.edges,
                    )
                end
                simple_to_nauty[simple_key] = nauty_key
            end
        end

        if first_simple_split !== nothing
            println("deduplicated witness split at order $order: $first_simple_split")
        end
        if first_simple_collapse !== nothing
            println("deduplicated witness collapse at order $order: $first_simple_collapse")
        end
        flush(stdout)

        @test length(nauty_to_multigraph) == expected
        @test length(multigraph_to_nauty) == expected
        println("rooted topology order $order -> $expected classes with exact multiplicities")
        flush(stdout)
    end
end
