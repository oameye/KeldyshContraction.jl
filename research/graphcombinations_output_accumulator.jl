using BenchmarkTools
using Test

import KeldyshContraction as KC

include(ENV["KC_CORPUS_DEFS"])

const OUTPUT_RESIDUAL_LABELS = Set([
    "boson/2-loop/none/none/sparse",
    "boson/2-loop/none/none/dense",
    "boson/2-loop/none/single/sparse",
    "boson/2-loop/none/pair/dense",
    "boson/2-loop/shell/none/dense",
    "boson/2-loop/shell/pair/dense",
    "boson/2-loop/pv/none/sparse",
    "fermion/2-loop/none/none/sparse",
    "fermion/2-loop/none/none/dense",
    "fermion/2-loop/shell/none/dense",
    "fermion/2-loop/shell/pair/dense",
    "fermion/2-loop/pv/none/sparse",
    "fermion/2-loop/pv/single/sparse",
])

function _workspace(expression)
    graph_capacity, loop_capacity = KC._loop_gc_capacities(expression)
    return KC._LoopGCQuotientWorkspace(graph_capacity, loop_capacity)
end

function prepared_current_output(
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

function _accumulate_kernel_atom!(
    out::Dict{KC.CollisionKernelSector{S},Dict{KC.OccupationMonomial{S},C}},
    sector::KC.CollisionKernelSector{S},
    monomial::KC.OccupationMonomial{S},
    coefficient::C,
) where {C<:Complex,S<:KC.Statistics}
    iszero(coefficient) && return out

    for (kinematic_monomial, kinematic_coefficient) in KC.kinematic_factor(sector)
        iszero(kinematic_coefficient) && continue
        unit_kinematic = KC.MomentumPolynomial(
            kinematic_monomial, one(KC.ComplexRationals)
        )
        unit_sector = KC.CollisionKernelSector{S}(
            KC.parameters(sector),
            KC.momentum_basis(sector),
            KC.external_wigner_momentum(sector),
            unit_kinematic,
            KC.frequency_support(sector),
        )
        KC._vanishes_under_loop_reflection(unit_sector, monomial) && continue

        scaled = convert(C, kinematic_coefficient) * coefficient
        iszero(scaled) && continue
        polynomial = if haskey(out, unit_sector)
            out[unit_sector]
        else
            created = Dict{KC.OccupationMonomial{S},C}()
            out[unit_sector] = created
            created
        end
        combined = get(polynomial, monomial, zero(C)) + scaled
        if iszero(combined)
            delete!(polynomial, monomial)
            isempty(polynomial) && delete!(out, unit_sector)
        else
            polynomial[monomial] = combined
        end
    end
    return out
end

function _materialize_accumulator(
    accumulator::Dict{
        KC.CollisionKernelSector{S},Dict{KC.OccupationMonomial{S},C}
    },
) where {C<:Number,S<:KC.Statistics}
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{C,S}}()
    sizehint!(out, length(accumulator))
    for (sector, polynomial) in accumulator
        isempty(polynomial) && continue
        terms = Pair{KC.OccupationMonomial{S},C}[]
        sizehint!(terms, length(polynomial))
        for (monomial, coefficient) in polynomial
            iszero(coefficient) || push!(terms, monomial => coefficient)
        end
        isempty(terms) && continue
        out[sector] = KC.OccupationPolynomial{C,S}(terms)
    end
    return out
end

function accumulated_output(
    expression::KC.OccupationReducedExpression{C,S,O,G,Ctx},
    workspace::KC._LoopGCQuotientWorkspace,
) where {C<:Number,S<:KC.Statistics,O,G,Ctx<:KC.AbstractWignerContext}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    accumulator = Dict{
        KC.CollisionKernelSector{S},Dict{KC.OccupationMonomial{S},D}
    }()

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
                _accumulate_kernel_atom!(
                    accumulator,
                    transformed_sector,
                    transformed_monomial,
                    transformed_coefficient,
                )
            end
        end
    end

    return KC.LoopQuotientedExpression{D,S,O,G,Ctx}(
        _materialize_accumulator(accumulator),
        KC.target_family(expression),
        KC.parameters(expression),
        KC.wigner_context(expression),
    )
end

cases = [
    label => expression for
    (label, expression) in corpus_cases() if label in OUTPUT_RESIDUAL_LABELS
]
workspaces = [_workspace(expression) for (_, expression) in cases]
println("direct output-accumulator residual corpus: ", length(cases), " cases")

@testset "direct output-accumulator semantics" begin
    @test length(cases) == length(OUTPUT_RESIDUAL_LABELS)
    for (index, (_, expression)) in enumerate(cases)
        workspace = workspaces[index]
        @test KC.loop_quotient_terms(accumulated_output(expression, workspace)) ==
            KC.loop_quotient_terms(prepared_current_output(expression, workspace))
    end
end

time_regressions = Tuple{String,Float64}[]
memory_regressions = Tuple{String,Int,Int}[]
current_ratios = Tuple{String,Float64}[]

@testset "direct output-accumulator residual performance" begin
    for (index, (label, expression)) in enumerate(cases)
        workspace = workspaces[index]
        nauty_quotient_loop_momenta(expression)
        prepared_current_output(expression, workspace)
        accumulated_output(expression, workspace)

        nauty_trial = @benchmark nauty_quotient_loop_momenta($expression) samples = 41 evals = 1
        current_trial = @benchmark prepared_current_output($expression, $workspace) samples = 41 evals = 1
        accumulator_trial = @benchmark accumulated_output($expression, $workspace) samples = 41 evals = 1
        nauty = median(nauty_trial)
        current = median(current_trial)
        accumulator = median(accumulator_trial)
        nauty_ratio = accumulator.time / nauty.time
        current_ratio = accumulator.time / current.time
        push!(current_ratios, (label, current_ratio))
        accumulator.time > nauty.time && push!(time_regressions, (label, nauty_ratio))
        accumulator.memory > nauty.memory &&
            push!(memory_regressions, (label, accumulator.memory, nauty.memory))
        println(
            label,
            ": accumulator=",
            round(accumulator.time / 1.0e3; digits=2),
            " μs/",
            accumulator.memory,
            " B; current=",
            round(current.time / 1.0e3; digits=2),
            " μs/",
            current.memory,
            " B; Nauty=",
            round(nauty.time / 1.0e3; digits=2),
            " μs/",
            nauty.memory,
            " B; ratios=",
            round(current_ratio; digits=3),
            "x current, ",
            round(nauty_ratio; digits=3),
            "x Nauty",
        )
    end

    println("confirmed time regressions: ", length(time_regressions))
    println("memory regressions: ", length(memory_regressions))
    println("accumulator/current ratios:")
    for (label, ratio) in sort(current_ratios; by=last, rev=true)
        println("  ", label, ": ", round(ratio; digits=3), "x")
    end
    if !isempty(time_regressions)
        println("remaining accumulator/Nauty time regressions:")
        for (label, ratio) in sort(time_regressions; by=last, rev=true)
            println("  ", label, ": ", round(ratio; digits=3), "x")
        end
    end
    if !isempty(memory_regressions)
        println("remaining accumulator/Nauty memory regressions:")
        for (label, accumulator_memory, nauty_memory) in memory_regressions
            println("  ", label, ": ", accumulator_memory, " B vs ", nauty_memory, " B")
        end
    end

    @test isempty(time_regressions)
    @test isempty(memory_regressions)
end
