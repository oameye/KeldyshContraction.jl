using BenchmarkTools
using KeldyshContraction
using Test

import KeldyshContraction as KC

include(joinpath(@__DIR__, "..", "benchmarks", "collision_reduction.jl"))

function nauty_quotient_loop_momenta(
    expression::KC.OccupationReducedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:KC.Statistics,O,G,Ctx<:KC.AbstractWignerContext}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{D,S}}()

    for (sector, polynomial) in KC.occupation_reduced_terms(expression)
        for (occupation_monomial, occupation_coefficient) in polynomial
            for (kinematic_monomial, kinematic_coefficient) in KC.kinematic_factor(sector)
                atom_sector = KC._kinematic_atom_sector(sector, kinematic_monomial)
                transform = KC._projective_canonical_loop_transform(
                    atom_sector, occupation_monomial
                )
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

function gc_oracle_workspace(
    sector::KC.ReducedCollisionSector, monomial::KC.OccupationMonomial
)
    basis = KC.momentum_basis(sector)
    external_index = KC._external_basis_index_noalloc(
        basis, KC.external_wigner_momentum(sector)
    )
    nloops = length(basis) - 1
    graph_capacity =
        3 + 3 * nloops + KC._projective_support_vertex_count(sector, external_index)

    for atom in monomial
        graph_capacity += 1 + KC._loop_incidence_vertex_count(atom.momentum, external_index)
    end

    kinematic_capacity = 0
    for (kinematic_monomial, _) in KC.kinematic_factor(sector)
        vertices = 1
        for component in kinematic_monomial
            vertices += KC._projective_component_vertex_count(component)
        end
        kinematic_capacity = max(kinematic_capacity, vertices)
    end
    graph_capacity += kinematic_capacity
    return KC._LoopGCQuotientWorkspace(graph_capacity, nloops)
end

function gc_requotient_terms(
    expression::KC.LoopQuotientedExpression{C,S}
) where {C<:Number,S<:KC.Statistics}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{D,S}}()

    for (sector, polynomial) in KC.loop_quotient_terms(expression)
        source_sector = reduced_sector(sector)
        for (monomial, coefficient) in polynomial
            workspace = gc_oracle_workspace(source_sector, monomial)
            transform = KC._graphcombinations_projective_canonical_loop_transform(
                source_sector, monomial, workspace
            )
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

function report_quotient_benchmark(nloops::Int, expression)
    nauty_quotient_loop_momenta(expression)
    KC.quotient_loop_momenta(expression)
    KC._loop_gc_capacities(expression)

    nauty_trial = @benchmark nauty_quotient_loop_momenta($expression) samples = 11 evals = 1
    gc_trial = @benchmark KC.quotient_loop_momenta($expression) samples = 11 evals = 1
    capacity_trial = @benchmark KC._loop_gc_capacities($expression) samples = 51 evals = 1
    nauty = median(nauty_trial)
    gc = median(gc_trial)
    capacity = median(capacity_trial)
    println(
        "$nloops-loop full quotient: Nauty ",
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
        " allocs; time ratio=",
        round(gc.time / nauty.time; digits=3),
        "; memory ratio=",
        round(gc.memory / nauty.memory; digits=3),
    )
    graph_capacity, loop_capacity = KC._loop_gc_capacities(expression)
    println(
        "$nloops-loop capacity scan: ",
        round(capacity.time / 1.0e3; digits=2),
        " μs / ",
        capacity.memory,
        " B / ",
        capacity.allocs,
        " allocs; graph_capacity=",
        graph_capacity,
        "; loop_capacity=",
        loop_capacity,
    )
    return nothing
end

@testset "production GC loop quotient vs retained Nauty oracle" begin
    for nloops in (2, 4)
        expression = benchmark_loop_quotient_fixture(nloops)
        nauty = nauty_quotient_loop_momenta(expression)
        gc = KC.quotient_loop_momenta(expression)
        gc_terms = KC.loop_quotient_terms(gc)

        @test gc_requotient_terms(gc) == gc_terms
        @test gc_requotient_terms(nauty) == gc_terms
        report_quotient_benchmark(nloops, expression)
    end
end
