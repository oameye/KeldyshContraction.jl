using BenchmarkTools
using KeldyshContraction
using Test

import GraphCombinations as GC
import KeldyshContraction as KC

const GC_SHA = "d7823f60d7f8a50b346ef68a0a0f927aedabec0c"
const OccupationPolynomial = KC.OccupationPolynomial
const OccupationAtom = KC.OccupationAtom

include(joinpath(@__DIR__, "..", "benchmarks", "collision_reduction.jl"))

function gc_result(edges::Vector{Tuple{Int,Int}}, colors::Vector{Int})
    graph = GC.DirectedGCGraph(
        Pair{Int,Int}[source => target for (source, target) in edges], length(colors)
    )
    return GC.canonicalize_directed(graph, colors)
end

function gc_loop_graph_data(
    sector::KC.ReducedCollisionSector{S}, monomial::KC.OccupationMonomial{S}
) where {S<:KC.Statistics}
    basis = KC.momentum_basis(sector)
    external = KC.external_wigner_momentum(sector)
    external_index = KC._external_basis_index(basis, external)
    loop_indices = KC._loop_basis_indices(basis, external)
    nloops = length(loop_indices)

    builder = KC._LoopCanonicalGraphBuilder()
    root = KC._add_loop_vertex!(builder, KC._loop_graph_color(1))
    pair_vertices = Vector{Int}(undef, nloops)
    positive_vertices = Vector{Int}(undef, nloops)
    negative_vertices = Vector{Int}(undef, nloops)
    loop_incidences = [Tuple{Int,KC.MomentumCoefficient}[] for _ in 1:nloops]

    for slot in 1:nloops
        pair = KC._add_loop_vertex!(builder, KC._loop_graph_color(2))
        positive = KC._add_loop_vertex!(builder, KC._loop_graph_color(3))
        negative = KC._add_loop_vertex!(builder, KC._loop_graph_color(3))
        pair_vertices[slot] = pair
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
    result = gc_result(builder.edges, labels)
    return result, pair_vertices, loop_incidences, basis, external, nloops
end

function gc_canonical_loop_transform(
    sector::KC.ReducedCollisionSector{S}, monomial::KC.OccupationMonomial{S}
) where {S<:KC.Statistics}
    result, pair_vertices, loop_incidences, basis, external, nloops = gc_loop_graph_data(
        sector, monomial
    )

    rank = GC.vertex_mapping(GC.canonical_relabeling(result))
    ordered_slots = sortperm(1:nloops; by=slot -> rank[pair_vertices[slot]])
    loop_permutation = zeros(Int, nloops)
    loop_signs = ones(Int, nloops)
    for (new_slot, old_slot) in enumerate(ordered_slots)
        loop_permutation[old_slot] = new_slot
        incidences = loop_incidences[old_slot]
        isempty(incidences) && continue
        canonical_first = first(sort!(copy(incidences); by=record -> rank[first(record)]))
        loop_signs[old_slot] = last(canonical_first) > 0 ? 1 : -1
    end
    return KC.loop_permutation_transform(
        basis, external, loop_permutation, loop_signs, zeros(Int, nloops)
    )
end

function gc_quotient_loop_momenta(
    expression::KC.OccupationReducedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:KC.Statistics,O,G,Ctx<:KC.AbstractWignerContext}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{D,S}}()

    for (sector, polynomial) in KC.occupation_reduced_terms(expression)
        for (occupation_monomial, occupation_coefficient) in polynomial
            for (kinematic_monomial, kinematic_coefficient) in KC.kinematic_factor(sector)
                atom_sector = KC._kinematic_atom_sector(sector, kinematic_monomial)
                transform = gc_canonical_loop_transform(atom_sector, occupation_monomial)
                transformed_sector, support_factor = KC._transform_kernel_sector(
                    atom_sector, transform
                )
                transformed_monomial = KC.transform_loop_momenta(
                    occupation_monomial, transform
                )
                transformed_coefficient =
                    convert(D, occupation_coefficient) *
                    convert(D, kinematic_coefficient) *
                    convert(D, support_factor)
                contribution = Pair{KC.OccupationMonomial{S},D}[transformed_monomial => transformed_coefficient]
                KC._push_kernel_polynomial!(
                    out, transformed_sector, KC.OccupationPolynomial{D,S}(contribution)
                )
            end
        end
    end

    return KC.LoopQuotientedExpression{D,S,O,G,Ctx}(
        out,
        KC.target_family(expression),
        KC.parameters(expression),
        KC.wigner_context(expression),
    )
end

function reduced_sector(sector::KC.CollisionKernelSector{S}) where {S<:KC.Statistics}
    return KC.ReducedCollisionSector{S}(
        KC.parameters(sector),
        KC.momentum_basis(sector),
        KC.external_wigner_momentum(sector),
        KC.kinematic_factor(sector),
        KC.frequency_support(sector),
    )
end

function gc_requotient_terms(
    expression::KC.LoopQuotientedExpression{C,S}
) where {C<:Number,S<:KC.Statistics}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{D,S}}()

    for (sector, polynomial) in KC.loop_quotient_terms(expression)
        source_sector = reduced_sector(sector)
        for (monomial, coefficient) in polynomial
            transform = gc_canonical_loop_transform(source_sector, monomial)
            transformed_sector, support_factor = KC._transform_kernel_sector(
                source_sector, transform
            )
            transformed_monomial = KC.transform_loop_momenta(monomial, transform)
            transformed_coefficient = convert(D, coefficient) * convert(D, support_factor)
            contribution = Pair{KC.OccupationMonomial{S},D}[transformed_monomial => transformed_coefficient]
            KC._push_kernel_polynomial!(
                out, transformed_sector, KC.OccupationPolynomial{D,S}(contribution)
            )
        end
    end
    return out
end

function report_equivalence(nloops::Int)
    expression = benchmark_loop_quotient_fixture(nloops)
    nauty = KC.quotient_loop_momenta(expression)
    gc = gc_quotient_loop_momenta(expression)
    nauty_terms = KC.loop_quotient_terms(nauty)
    gc_terms = KC.loop_quotient_terms(gc)
    normalized_nauty = gc_requotient_terms(nauty)
    normalized_gc = gc_requotient_terms(gc)

    println("GC loop equivalence probe ($GC_SHA): $nloops loops")
    println("direct exact equality: ", gc_terms == nauty_terms)
    println("Nauty sectors: ", length(nauty_terms), "; GC sectors: ", length(gc_terms))
    println("GC re-quotient idempotent: ", normalized_gc == gc_terms)
    println("orbit-normalized equality: ", normalized_nauty == gc_terms)

    @test normalized_gc == gc_terms
    @test normalized_nauty == gc_terms
    return expression
end

function report_benchmark(nloops::Int, expression)
    KC.quotient_loop_momenta(expression)
    gc_quotient_loop_momenta(expression)
    nauty_trial = @benchmark KC.quotient_loop_momenta($expression) samples = 7 evals = 1
    gc_trial = @benchmark gc_quotient_loop_momenta($expression) samples = 7 evals = 1
    nauty = median(nauty_trial)
    gc = median(gc_trial)
    println(
        "$nloops-loop warmed median: Nauty ",
        round(nauty.time / 1.0e3; digits=2),
        " μs / ",
        nauty.memory,
        " B / ",
        nauty.allocs,
        " allocs; GC ",
        round(gc.time / 1.0e3; digits=2),
        " μs / ",
        gc.memory,
        " B / ",
        gc.allocs,
        " allocs",
    )
    return nothing
end

@testset "GC/Nauty loop-orbit equivalence" begin
    for nloops in (2, 4)
        expression = report_equivalence(nloops)
        report_benchmark(nloops, expression)
    end
end
