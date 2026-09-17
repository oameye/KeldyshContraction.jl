using BenchmarkTools
using Test

import KeldyshContraction as KC

include(ENV["KC_CORPUS_DEFS"])

function reusable_gc_quotient_loop_momenta(
    expression::KC.OccupationReducedExpression{C,S,O,G,Ctx},
    workspace::KC._LoopGCQuotientWorkspace,
) where {C<:Number,S<:KC.Statistics,O,G,Ctx<:KC.AbstractWignerContext}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{D,S}}()

    for (sector, polynomial) in KC.occupation_reduced_terms(expression)
        for (occupation_monomial, occupation_coefficient) in polynomial
            for (kinematic_monomial, kinematic_coefficient) in KC.kinematic_factor(sector)
                atom_sector = KC._kinematic_atom_sector(sector, kinematic_monomial)
                transform = KC._graphcombinations_projective_canonical_loop_transform(
                    atom_sector, occupation_monomial, workspace
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
                contribution = Pair{KC.OccupationMonomial{S},D}[
                    transformed_monomial => transformed_coefficient
                ]
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

function _prepared_workspace(expression)
    graph_capacity, loop_capacity = KC._loop_gc_capacities(expression)
    return KC._LoopGCQuotientWorkspace(graph_capacity, loop_capacity)
end

function _measure_prepared_pair_15(expression, workspace)
    nauty_quotient_loop_momenta(expression)
    reusable_gc_quotient_loop_momenta(expression, workspace)
    nauty_trial = @benchmark nauty_quotient_loop_momenta($expression) samples = 15 evals = 1
    gc_trial = @benchmark reusable_gc_quotient_loop_momenta(
        $expression, $workspace
    ) samples = 15 evals = 1
    return median(nauty_trial), median(gc_trial)
end

function _measure_prepared_pair_41(expression, workspace)
    nauty_trial = @benchmark nauty_quotient_loop_momenta($expression) samples = 41 evals = 1
    gc_trial = @benchmark reusable_gc_quotient_loop_momenta(
        $expression, $workspace
    ) samples = 41 evals = 1
    return median(nauty_trial), median(gc_trial)
end

function _confirmed_prepared_measurement(expression, workspace)
    nauty, gc = _measure_prepared_pair_15(expression, workspace)
    if gc.time > nauty.time || gc.memory > nauty.memory
        nauty, gc = _measure_prepared_pair_41(expression, workspace)
    end
    return nauty, gc
end

cases = corpus_cases()
workspaces = [_prepared_workspace(expression) for (_, expression) in cases]
println(
    "prepared packed-GC/Nauty corpus: ",
    length(cases),
    " generated production-shaped cases",
)

semantic_failures = String[]
time_regressions = Tuple{String,Float64}[]
memory_regressions = Tuple{String,Int,Int}[]
worst_ratio = ("", 0.0)
best_ratio = ("", Inf)

@testset "prepared packed-GC semantics" begin
    for (index, (label, expression)) in enumerate(cases)
        workspace = workspaces[index]
        prepared = reusable_gc_quotient_loop_momenta(expression, workspace)
        current = KC.quotient_loop_momenta(expression)
        prepared_terms = KC.loop_quotient_terms(prepared)
        exact =
            prepared_terms == KC.loop_quotient_terms(current) &&
            gc_requotient_terms(nauty_quotient_loop_momenta(expression)) == prepared_terms
        if !exact
            push!(semantic_failures, label)
        end
        @test exact
    end
end

@testset "prepared packed-GC performance" begin
    global worst_ratio, best_ratio
    for (index, (label, expression)) in enumerate(cases)
        workspace = workspaces[index]
        nauty, gc = _confirmed_prepared_measurement(expression, workspace)
        ratio = gc.time / nauty.time
        if ratio > worst_ratio[2]
            worst_ratio = (label, ratio)
        end
        if ratio < best_ratio[2]
            best_ratio = (label, ratio)
        end
        if gc.time > nauty.time
            push!(time_regressions, (label, ratio))
        end
        if gc.memory > nauty.memory
            push!(memory_regressions, (label, gc.memory, nauty.memory))
        end
        println(
            lpad(index, 3),
            "/",
            length(cases),
            " ",
            label,
            ": time=",
            round(gc.time / 1.0e3; digits=2),
            "/",
            round(nauty.time / 1.0e3; digits=2),
            " μs prepared-GC/Nauty (",
            round(ratio; digits=3),
            "x), memory=",
            gc.memory,
            "/",
            nauty.memory,
            " B",
        )
    end

    println("semantic failures: ", length(semantic_failures))
    println("confirmed time regressions: ", length(time_regressions))
    println("memory regressions: ", length(memory_regressions))
    println("best time ratio: ", round(best_ratio[2]; digits=3), "x at ", best_ratio[1])
    println("worst time ratio: ", round(worst_ratio[2]; digits=3), "x at ", worst_ratio[1])

    if !isempty(time_regressions)
        println("confirmed prepared-GC time regressions:")
        for (label, ratio) in sort(time_regressions; by=last, rev=true)
            println("  ", label, ": ", round(ratio; digits=3), "x")
        end
    end
    if !isempty(memory_regressions)
        println("prepared-GC memory regressions:")
        for (label, gc_memory, nauty_memory) in memory_regressions
            println("  ", label, ": ", gc_memory, " B vs ", nauty_memory, " B")
        end
    end

    @test isempty(time_regressions)
    @test isempty(memory_regressions)
end
