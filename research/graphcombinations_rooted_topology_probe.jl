using KeldyshContraction
using Test

import GraphCombinations as GC
import KeldyshContraction as KC

const GC_SHA = "e0e5d7dc4e14f4fc356ffeee3809fdb5b47fcf12"

function gc_legacy_graph(vs)
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
        vertices in edges || push!(edges, vertices)
    end

    return GC.DirectedGCGraph(
        Pair{Int,Int}[source => target for (source, target) in edges], max_label
    )
end

function certify_gc_graph_kernel(vs)
    graph = gc_legacy_graph(vs)
    reference = GC.canonicalize_directed(graph)

    buffer = GC.DirectedCanonicalizationBuffer(graph.num_vertices)
    workspace = GC.DirectedCanonicalizationWorkspace(graph.num_vertices)
    GC.canonicalize_directed!(buffer, workspace, graph)

    @test Tuple(buffer.canonical_multiplicities) ==
        Tuple(GC.canonical_graph(reference).multiplicities)
    @test buffer.old_to_canonical == GC.vertex_mapping(GC.canonical_relabeling(reference))
    return nothing
end

@qfields rooted_probe_ϕ::Boson
c, q = rooted_probe_ϕ[Classical], rooted_probe_ϕ[Quantum]

println("GC rooted topology compatibility probe ($GC_SHA)")
flush(stdout)

@testset "legacy Nauty topology contract and GC graph exactness" begin
    elastic = -(
        1 // 2 * (c^2 + q^2) * bar(c) * bar(q) + 1 // 2 * c * q * (bar(c)^2 + bar(q)^2)
    )
    L = InteractionLagrangian(elastic)
    cases = ((1, 3, 1), (2, 5, 3), (3, 7, 11), (4, 9, 59))

    for (order, edge_count, expected) in cases
        G = DressedPropagator(L, Val(order), Val(edge_count))
        component = KC.topologies(G.keldysh)

        # The 1/3/11/59 vectors are a historical KC compatibility contract derived
        # from Nauty's canonical-label witness. GC deliberately does not reproduce
        # another backend's incidental vertex numbering; issue GC #151 defines that
        # numbering as consumer metadata rather than a generic GC correctness oracle.
        @test length(keys(component)) == expected

        checked = 0
        for diagrams in values(component)
            for diagram in diagrams
                contractions = KC.Contraction{Boson}[
                    (edge.out, edge.in) for edge in KC.contractions(diagram)
                ]
                certify_gc_graph_kernel(contractions)
                checked += 1
            end
        end

        @test checked > 0
        println(
            "topology order $order: legacy Nauty=$expected classes; ",
            "GC reference/workspace exact on $checked routed graphs",
        )
        flush(stdout)
    end
end
