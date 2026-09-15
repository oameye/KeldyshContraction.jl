using BenchmarkTools
using KeldyshContraction
using Test

import GraphCombinations as GC
import KeldyshContraction as KC

const GC_SHA = "a09d6fdbd86bb4accae6ffbd3be1e0c801c49c05"

include(joinpath(@__DIR__, "..", "benchmarks", "collision_reduction.jl"))

function gc_loop_graph(
    sector::KC.ReducedCollisionSector{S}, monomial::KC.OccupationMonomial{S}
) where {S<:KC.Statistics}
    basis = KC.momentum_basis(sector)
    external = KC.external_wigner_momentum(sector)
    external_index = KC._external_basis_index(basis, external)
    loop_indices = KC._loop_basis_indices(basis, external)
    nloops = length(loop_indices)

    builder = KC._LoopCanonicalGraphBuilder()
    root = KC._add_loop_vertex!(builder, KC._loop_graph_color(1))
    positive_vertices = Vector{Int}(undef, nloops)
    negative_vertices = Vector{Int}(undef, nloops)
    loop_incidences = [Tuple{Int,KC.MomentumCoefficient}[] for _ in 1:nloops]

    for slot in 1:nloops
        pair = KC._add_loop_vertex!(builder, KC._loop_graph_color(2))
        positive = KC._add_loop_vertex!(builder, KC._loop_graph_color(3))
        negative = KC._add_loop_vertex!(builder, KC._loop_graph_color(3))
        positive_vertices[slot] = positive
        negative_vertices[slot] = negative
        KC._add_loop_edge!(builder, root, pair)
        KC._add_loop_edge!(builder, pair, positive)
        KC._add_loop_edge!(builder, pair, negative)
    end

    KC._add_loop_semantics!(
        builder,
        root,
        sector,
        monomial,
        external_index,
        loop_indices,
        positive_vertices,
        negative_vertices,
        loop_incidences,
    )

    color_classes = sort!(unique(copy(builder.colors)))
    labels = Int[searchsortedfirst(color_classes, color) for color in builder.colors]
    edges = Pair{Int,Int}[source => target for (source, target) in builder.edges]
    graph = GC.DirectedGCGraph(edges, length(labels))
    return graph, labels, builder
end

function gc_loop_graphs(expression)
    graphs = Tuple{GC.DirectedGCGraph,Vector{Int},KC._LoopCanonicalGraphBuilder}[]
    for (sector, polynomial) in KC.occupation_reduced_terms(expression)
        for (occupation_monomial, _) in polynomial
            for (kinematic_monomial, _) in KC.kinematic_factor(sector)
                atom_sector = KC._kinematic_atom_sector(sector, kinematic_monomial)
                push!(graphs, gc_loop_graph(atom_sector, occupation_monomial))
            end
        end
    end
    return graphs
end

function certify_workspace(graph::GC.DirectedGCGraph, labels::Vector{Int})
    expected = GC.canonicalize_directed(graph, labels)
    workspace = GC.DirectedCanonicalizationWorkspace(graph.num_vertices)
    buffer = GC.DirectedCanonicalizationBuffer(graph.num_vertices)

    @test @inferred(GC.canonicalize_directed!(buffer, workspace, graph, labels)) === buffer
    @test GC.canonical_graph(buffer) == GC.canonical_graph(expected)
    @test GC.canonical_automorphism_order(buffer) == GC.canonical_automorphism_order(expected)

    mapping = GC.vertex_mapping(GC.canonical_relabeling(expected))
    @test all(GC.canonical_rank(buffer, vertex) == mapping[vertex] for vertex in eachindex(mapping))
    @test all(
        GC.original_vertex(buffer, GC.canonical_rank(buffer, vertex)) == vertex for
        vertex in eachindex(mapping)
    )

    GC.canonicalize_directed!(buffer, workspace, graph, labels)
    @test @allocated(GC.canonicalize_directed!(buffer, workspace, graph, labels)) == 0
    return workspace, buffer
end

function report_kernel_benchmark(nloops::Int, graph, labels, builder, workspace, buffer)
    nauty_graph = KC._loop_graph(builder)
    KC.NautyGraphs.canonical_permutation(nauty_graph)
    GC.canonicalize_directed!(buffer, workspace, graph, labels)

    nauty_trial = @benchmark KC.NautyGraphs.canonical_permutation($nauty_graph) samples = 9 evals = 1
    gc_trial = @benchmark GC.canonicalize_directed!(
        $buffer, $workspace, $graph, $labels
    ) samples = 9 evals = 1
    nauty = median(nauty_trial)
    gc = median(gc_trial)
    println(
        "$nloops-loop canonical-label kernel ($GC_SHA): Nauty ",
        round(nauty.time / 1.0e3; digits=2),
        " μs / ",
        nauty.memory,
        " B / ",
        nauty.allocs,
        " allocs; GC workspace ",
        round(gc.time / 1.0e3; digits=2),
        " μs / ",
        gc.memory,
        " B / ",
        gc.allocs,
        " allocs",
    )
    return nothing
end

@testset "GC production workspace on KC loop graphs" begin
    for nloops in (2, 4)
        expression = benchmark_loop_quotient_fixture(nloops)
        graphs = gc_loop_graphs(expression)
        @test !isempty(graphs)

        for (graph, labels, _) in graphs
            certify_workspace(graph, labels)
        end

        graph, labels, builder = first(graphs)
        workspace, buffer = certify_workspace(graph, labels)
        report_kernel_benchmark(nloops, graph, labels, builder, workspace, buffer)
    end
end
