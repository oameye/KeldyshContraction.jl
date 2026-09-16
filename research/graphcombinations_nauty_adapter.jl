using BenchmarkTools
using KeldyshContraction
using Test

import GraphCombinations as GC
import KeldyshContraction as KC

const GC_SHA = "e0e5d7dc4e14f4fc356ffeee3809fdb5b47fcf12"
const Out = KC.Out
const In = KC.In
const Bulk = KC.Bulk

function gc_result(edges::Vector{Tuple{Int,Int}}, colors::Vector{Int})
    graph = GC.DirectedGCGraph(
        Pair{Int,Int}[source => target for (source, target) in edges], length(colors)
    )
    result = GC.canonicalize_directed(graph, colors)

    workspace = GC.DirectedCanonicalizationWorkspace(graph.num_vertices)
    buffer = GC.DirectedCanonicalizationBuffer(graph.num_vertices)
    GC.canonicalize_directed!(buffer, workspace, graph, colors)
    @test GC.canonical_graph(buffer) == GC.canonical_graph(result)
    @test GC.canonical_automorphism_order(buffer) == GC.canonical_automorphism_order(result)
    @test buffer.old_to_canonical == GC.vertex_mapping(GC.canonical_relabeling(result))
    return result
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

function gc_topology_result(vs, graph_positions=KC.canonicalization_positions(vs))
    edges = unique(gc_position_edges(vs, graph_positions))
    colors = KC.position_labels(graph_positions)
    return gc_result(edges, colors)
end

function gc_physical_result(vs, graph_positions=KC.canonicalization_positions(vs))
    direct_edges = gc_position_edges(vs, graph_positions)
    simple = length(unique(direct_edges)) == length(direct_edges)
    if simple && KC.uniform_coloring(vs)
        return gc_result(direct_edges, KC.position_labels(graph_positions))
    end

    physical_colors = KC.propagator_colors(vs)
    npositions = length(graph_positions)
    labels = Vector{Int}(undef, npositions + length(vs))
    copyto!(labels, 1, KC.position_labels(graph_positions), 1, npositions)
    edges = Tuple{Int,Int}[]
    sizehint!(edges, 2 * length(vs))
    for (i, item) in enumerate(vs)
        color_index = searchsortedfirst(physical_colors, KC.propagator_color(item))
        edge_vertex = npositions + i
        labels[edge_vertex] = 3 + color_index
        out, in = KC.positions(item)
        push!(edges, (KC.position_vertex(graph_positions, out), edge_vertex))
        push!(edges, (edge_vertex, KC.position_vertex(graph_positions, in)))
    end
    return gc_result(edges, labels)
end

function gc_canonicalize(vs::Vector{T}) where {T}
    isempty(vs) && return copy(vs)
    graph_positions = KC.canonicalization_positions(vs)
    result = gc_physical_result(vs, graph_positions)
    old_to_canonical = GC.vertex_mapping(GC.canonical_relabeling(result))
    canonical_to_old = sortperm(old_to_canonical)
    mapping = KC.make_permutation_dict(canonical_to_old, graph_positions, vs)
    return T[KC.relabel_bulk_positions(item, mapping) for item in vs]
end

@qfields gc_adapter_ϕ::Boson
@qfields gc_adapter_χ::Boson

c, q = gc_adapter_ϕ[Classical], gc_adapter_ϕ[Quantum]
χc, χq = gc_adapter_χ[Classical], gc_adapter_χ[Quantum]

function as_contractions(vs)
    return KC.Contraction{Boson}[KC.Contraction(item) for item in vs]
end

println("GC physical canonicalization oracle ($GC_SHA)")
flush(stdout)

@testset "GraphCombinations physical KC oracle" begin
    ring = as_contractions([
        (c(Out()), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(Bulk(2))),
        (c(Bulk(2)), bar(q)(Bulk(3))),
        (c(Bulk(3)), bar(q)(Bulk(4))),
        (c(Bulk(4)), bar(q)(Bulk(1))),
        (c(Bulk(4)), bar(q)(In())),
    ])
    self_loop = as_contractions([
        (c(Out()), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(In())),
    ])
    repeated = as_contractions([
        (c(Out()), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(Bulk(2))),
        (c(Bulk(1)), bar(q)(Bulk(2))),
        (c(Bulk(2)), bar(q)(In())),
    ])
    colored = as_contractions([
        (c(Out()), bar(q)(Bulk(1))),
        (χc(Bulk(1)), bar(χq)(Bulk(2))),
        (c(Bulk(2)), bar(q)(In())),
    ])
    symmetric = as_contractions([
        (c(Bulk(1)), bar(q)(Bulk(2))),
        (c(Bulk(2)), bar(q)(Bulk(3))),
        (c(Bulk(3)), bar(q)(Bulk(1))),
    ])

    for fixture in (ring, self_loop, repeated, colored, symmetric)
        gc_fixture = gc_canonicalize(fixture)
        nauty_fixture = KC.canonicalize(fixture)
        @test gc_canonicalize(nauty_fixture) == gc_fixture
        @test KC.canonicalize(gc_fixture) == nauty_fixture

        graph_positions = KC.canonicalization_positions(fixture)
        gc_topology = gc_topology_result(fixture, graph_positions)
        _, _, _, nauty_automorphisms = KC.canonicalization_permutations(
            fixture, graph_positions
        )
        @test GC.canonical_automorphism_order(gc_topology) == nauty_automorphisms.n
    end
    @test GC.canonical_automorphism_order(gc_topology_result(symmetric)) == 3

    KC.canonicalize(ring)
    gc_canonicalize(ring)
    nauty_trial = @benchmark KC.canonicalize($ring) samples = 7 evals = 1
    gc_trial = @benchmark gc_canonicalize($ring) samples = 7 evals = 1
    nauty = median(nauty_trial)
    gc = median(gc_trial)
    println(
        "propagator ring warmed median: Nauty ",
        round(nauty.time / 1.0e3; digits=2),
        " μs / ",
        nauty.memory,
        " B / ",
        nauty.allocs,
        " allocs; GC adapter ",
        round(gc.time / 1.0e3; digits=2),
        " μs / ",
        gc.memory,
        " B / ",
        gc.allocs,
        " allocs",
    )
end
