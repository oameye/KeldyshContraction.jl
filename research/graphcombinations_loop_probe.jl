using BenchmarkTools
using KeldyshContraction

import GraphCombinations as GC
import KeldyshContraction as KC

const GC_SHA = "ec5f84d4b1b1f21b2dacdb00ce4ed9af6bd68fa8"
const OccupationPolynomial = KC.OccupationPolynomial
const OccupationAtom = KC.OccupationAtom

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
    graph = GC.DirectedGCGraph(
        Pair{Int,Int}[source => target for (source, target) in builder.edges],
        length(labels),
    )
    return graph, labels, length(builder.edges)
end

function color_cell_sizes(labels)
    counts = Dict{Int,Int}()
    for label in labels
        counts[label] = get(counts, label, 0) + 1
    end
    return sort!(collect(values(counts)); rev=true)
end

function relabeling_group_order(cells)
    order = big(1)
    for count in cells
        order *= factorial(big(count))
    end
    return order
end

function loop_records(expression)
    records = NamedTuple[]
    for (sector, polynomial) in KC.occupation_reduced_terms(expression)
        for (monomial, _) in polynomial
            for (kinematic_monomial, _) in KC.kinematic_factor(sector)
                atom_sector = KC._kinematic_atom_sector(sector, kinematic_monomial)
                graph, labels, edge_count = gc_loop_graph(atom_sector, monomial)
                cells = color_cell_sizes(labels)
                push!(
                    records,
                    (
                        group_order=relabeling_group_order(cells),
                        vertices=length(labels),
                        edges=edge_count,
                        cells=cells,
                        graph=graph,
                        labels=labels,
                        sector=atom_sector,
                        monomial=monomial,
                    ),
                )
            end
        end
    end
    return records
end

function report_inventory(nloops)
    expression = benchmark_loop_quotient_fixture(nloops)
    records = loop_records(expression)
    isempty(records) && error("loop fixture produced no canonicalization records")
    worst_index = argmax(record.group_order for record in records)
    worst = records[worst_index]
    println("GC loop probe ($GC_SHA): $nloops loops")
    println("integrand atoms: ", length(records))
    println("worst semantic graph vertices/edges: ", worst.vertices, "/", worst.edges)
    println("worst color-cell sizes: ", worst.cells)
    println("exact enumerated relabeling group order: ", worst.group_order)
    flush(stdout)
    return worst
end

function report_measurement(nloops)
    worst = report_inventory(nloops)

    if nloops == 2
        KC._canonical_loop_transform(worst.sector, worst.monomial)
        GC.canonicalize_directed(worst.graph, worst.labels)
    end

    nauty = @timed KC._canonical_loop_transform(worst.sector, worst.monomial)
    println(
        "Nauty transform: ",
        round(nauty.time * 1.0e3; digits=3),
        " ms / ",
        nauty.bytes,
        " B",
    )
    flush(stdout)

    gc = @timed GC.canonicalize_directed(worst.graph, worst.labels)
    println(
        "GC canonicalize_directed: ",
        round(gc.time * 1.0e3; digits=3),
        " ms / ",
        gc.bytes,
        " B",
    )
    flush(stdout)
    return nothing
end

length(ARGS) == 2 ||
    error("usage: graphcombinations_loop_probe.jl <nloops> <inventory|measure>")
nloops = parse(Int, ARGS[1])
mode = ARGS[2]
if mode == "inventory"
    report_inventory(nloops)
elseif mode == "measure"
    report_measurement(nloops)
else
    error("unknown loop probe mode: $mode")
end
