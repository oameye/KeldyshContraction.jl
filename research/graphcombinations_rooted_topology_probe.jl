using KeldyshContraction
using Test

import GraphCombinations as GC
import KeldyshContraction as KC

const GC_SHA = "ec5f84d4b1b1f21b2dacdb00ce4ed9af6bd68fa8"
const Out = KC.Out
const In = KC.In
const Bulk = KC.Bulk

function gc_result(edges::Vector{Tuple{Int,Int}}, colors::Vector{Int})
    graph = GC.DirectedGCGraph(
        Pair{Int,Int}[source => target for (source, target) in edges], length(colors)
    )
    return GC.canonicalize_directed(graph, colors)
end

function gc_position_edges(vs, graph_positions)
    edges = Tuple{Int,Int}[]
    sizehint!(edges, length(vs))
    for item in vs
        out, in = KC.positions(item)
        push!(
            edges,
            (
                KC.position_vertex(graph_positions, out),
                KC.position_vertex(graph_positions, in),
            ),
        )
    end
    return edges
end

function gc_rooted_topology_mapping(old_to_canonical, graph_positions)
    canonical_to_old = sortperm(old_to_canonical)
    canonical_bulk = Int[
        original_vertex for
        original_vertex in canonical_to_old if KC.is_bulk(graph_positions[original_vertex])
    ]
    mapping = Dict{KC.Position,KC.Position}()
    for (bulk_index, original_vertex) in enumerate(canonical_bulk)
        mapping[graph_positions[original_vertex]] = Bulk(bulk_index)
    end
    return mapping
end

function gc_rooted_topology(vs, ::Val{E2}) where {E2}
    isempty(vs) && return KC.bulk_multiplicity(Tuple{Int8,Int8}[], Val(E2))

    graph_positions = KC.canonicalization_positions(vs)
    edges = unique(gc_position_edges(vs, graph_positions))
    colors = KC.position_labels(graph_positions)
    result = gc_result(edges, colors)

    old_to_canonical = GC.vertex_mapping(GC.canonical_relabeling(result))
    mapping = gc_rooted_topology_mapping(old_to_canonical, graph_positions)
    topology_edges = Tuple{Int8,Int8}[
        KC.integer_positions(KC.relabel_bulk_positions(item, mapping)) for item in vs
    ]
    return KC.bulk_multiplicity(topology_edges, Val(E2))
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

        nauty_to_gc = Dict{Any,Any}()
        gc_to_nauty = Dict{Any,Any}()
        for (key, diagrams) in component
            nauty_key = Tuple(key)
            rooted_keys = Set{Any}()
            for diagram in diagrams
                contractions = KC.Contraction{Boson}[
                    (edge.out, edge.in) for edge in KC.contractions(diagram)
                ]
                push!(
                    rooted_keys, Tuple(gc_rooted_topology(contractions, Val(length(key))))
                )
            end

            @test length(rooted_keys) == 1
            gc_key = only(rooted_keys)
            @test get!(nauty_to_gc, nauty_key, gc_key) == gc_key
            @test get!(gc_to_nauty, gc_key, nauty_key) == nauty_key
        end

        @test length(nauty_to_gc) == expected
        @test length(gc_to_nauty) == expected
        println("rooted topology order $order -> $expected classes")
        flush(stdout)
    end
end
