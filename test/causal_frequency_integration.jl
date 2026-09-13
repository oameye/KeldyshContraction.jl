using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields causal_expression_ϕ::Boson

function causal_expression_energies()
    basis = KC.MomentumBasis(3)
    momenta = ntuple(i -> KC.basis_momentum(basis, i), 3)
    return map(
        momentum -> EnergyForm(DispersionAtom(causal_expression_ϕ, momentum)), momenta
    )
end

function causal_expression_denominator(coefficients::NTuple{N,Int}, energy, infinitesimal) where {N}
    return KC.CausalFrequencyDenominator(
        [coefficient // 1 for coefficient in coefficients], energy, infinitesimal // 1
    )
end
function causal_expression_denominator(coefficient::Int, energy, infinitesimal)
    return causal_expression_denominator((coefficient,), energy, infinitesimal)
end

function one_term_expression(coefficient, denominators)
    term = KC.CausalFrequencyTerm(coefficient, denominators)
    return KC.CausalFrequencyExpression([term])
end

@testset "expression-level causal frequency integration" begin
    Eₐ, Eᵦ, Eᵧ = causal_expression_energies()

    @testset "global marginal-tail cancellation is contour safe" begin
        lower_a = causal_expression_denominator(1, -Eₐ, 1)
        lower_b = causal_expression_denominator(1, -Eᵦ, 1)
        term_a = KC.CausalFrequencyTerm(1 // 1, [lower_a])
        term_b = KC.CausalFrequencyTerm(-1 // 1, [lower_b])

        @test KC.frequency_decay_lower_bound(KC.CausalFrequencyExpression([term_a]), 1) == 1
        @test KC.frequency_decay_lower_bound(KC.CausalFrequencyExpression([term_b]), 1) == 1

        expression = KC.CausalFrequencyExpression([term_a, term_b])
        @test KC.causal_contour_safe_at_infinity(expression, 1)

        result = @inferred KC.integrate_causal_frequency_expression(expression, 1)
        @test KC.causal_frequency_integration_kind(result) === KC.CausalFrequencyIntegrated
        @test isempty(KC.causal_frequency_integration_expression(result))
    end

    @testset "global cancellation can have a nonzero residue sum" begin
        lower = causal_expression_denominator(1, -Eₐ, 1)
        upper = causal_expression_denominator(1, -Eᵦ, -1)
        expression = KC.CausalFrequencyExpression([
            KC.CausalFrequencyTerm(1 // 1, [lower]),
            KC.CausalFrequencyTerm(-1 // 1, [upper]),
        ])

        result = @inferred KC.integrate_causal_frequency_expression(expression, 1)
        @test KC.causal_frequency_integration_kind(result) === KC.CausalFrequencyIntegrated
        integrated = KC.causal_frequency_integration_expression(result)
        terms = KC.causal_frequency_terms(integrated)
        @test length(terms) == 1
        @test isempty(KC.causal_frequency_denominators(only(terms)))
        @test KC.causal_frequency_coefficient(only(terms)) == -im
    end

    @testset "same-half-plane regular expression vanishes" begin
        expression = one_term_expression(
            1 // 1,
            [
                causal_expression_denominator(1, -Eₐ, 1),
                causal_expression_denominator(1, -Eᵦ, 1),
            ],
        )
        result = @inferred KC.integrate_causal_frequency_expression(expression, 1)
        @test KC.causal_frequency_integration_kind(result) === KC.CausalFrequencyIntegrated
        @test isempty(KC.causal_frequency_integration_expression(result))
    end

    @testset "constant term is a typed infinity divergence" begin
        expression = one_term_expression(1 // 1, KC.CausalFrequencyDenominator{Boson}[])
        result = @inferred KC.integrate_causal_frequency_expression(expression, 1)
        @test KC.causal_frequency_integration_kind(result) ===
            KC.CausalFrequencyDivergentAtInfinity
    end

    @testset "marginal infinity is deferred without exception" begin
        expression = one_term_expression(1 // 1, [causal_expression_denominator(1, -Eₐ, 1)])
        result = @inferred KC.integrate_causal_frequency_expression(expression, 1)
        @test KC.causal_frequency_integration_kind(result) ===
            KC.CausalFrequencyMarginalAtInfinity
        @test KC.causal_frequency_integration_expression(result) ==
            KC._complex_causal_frequency_expression(expression)
    end

    @testset "zero prescription is deferred" begin
        expression = one_term_expression(
            1 // 1,
            [
                causal_expression_denominator(1, -Eₐ, 0),
                causal_expression_denominator(1, -Eᵦ, 1),
            ],
        )
        result = @inferred KC.integrate_causal_frequency_expression(expression, 1)
        @test KC.causal_frequency_integration_kind(result) ===
            KC.CausalFrequencyZeroPrescription
    end

    @testset "upper repeated pole is deferred" begin
        upper = causal_expression_denominator(1, -Eₐ, -1)
        expression = one_term_expression(1 // 1, [upper, upper])
        result = @inferred KC.integrate_causal_frequency_expression(expression, 1)
        @test KC.causal_frequency_integration_kind(result) ===
            KC.CausalFrequencyRepeatedPole
    end

    @testset "multi-frequency reduction is order invariant" begin
        lower_1 = causal_expression_denominator((1, 0), -Eₐ, 1)
        upper_1 = causal_expression_denominator((1, 0), -Eᵦ, -1)
        lower_2 = causal_expression_denominator((0, 1), -Eᵦ, 1)
        upper_2 = causal_expression_denominator((0, 1), -Eᵧ, -1)
        expression = KC.CausalFrequencyExpression([
            KC.CausalFrequencyTerm(1 // 1, [lower_1, lower_2]),
            KC.CausalFrequencyTerm(-1 // 1, [lower_1, upper_2]),
            KC.CausalFrequencyTerm(-1 // 1, [upper_1, lower_2]),
            KC.CausalFrequencyTerm(1 // 1, [upper_1, upper_2]),
        ])

        canonical_plan = @inferred KC.canonical_causal_frequency_reduction_plan(Val(2))
        reverse_plan = @inferred KC.CausalFrequencyReductionPlan((2, 1))
        @test KC.causal_frequency_reduction_order(canonical_plan) == (1, 2)
        @test KC.causal_frequency_reduction_order(reverse_plan) == (2, 1)
        @test_throws ArgumentError KC.CausalFrequencyReductionPlan((1, 1))

        forward = @inferred KC.reduce_causal_frequency_expression(expression, canonical_plan)
        reverse = @inferred KC.reduce_causal_frequency_expression(expression, reverse_plan)
        @test KC.causal_frequency_reduction_kind(forward) === KC.CausalFrequencyIntegrated
        @test KC.causal_frequency_reduction_kind(reverse) === KC.CausalFrequencyIntegrated
        @test KC.causal_frequency_reduction_steps(forward) == 2
        @test KC.causal_frequency_reduction_steps(reverse) == 2
        @test KC.causal_frequency_reduction_blocked_frequency(forward) == 0
        @test KC.causal_frequency_reduction_blocked_frequency(reverse) == 0
        @test KC.causal_frequency_reduction_expression(forward) ==
            KC.causal_frequency_reduction_expression(reverse)

        terms = KC.causal_frequency_terms(KC.causal_frequency_reduction_expression(forward))
        @test length(terms) == 1
        @test isempty(KC.causal_frequency_denominators(only(terms)))
        @test KC.causal_frequency_coefficient(only(terms)) == -1
    end

    @testset "coupled affine poles are order invariant" begin
        upper_1 = causal_expression_denominator((1, 0), -Eₐ, -1)
        upper_2 = causal_expression_denominator((0, 1), -Eᵦ, -1)
        coupled_lower = causal_expression_denominator((1, 1), -Eᵧ, 1)
        expression = one_term_expression(1 // 1, [upper_1, upper_2, coupled_lower])

        forward = @inferred KC.reduce_causal_frequency_expression(
            expression, KC.CausalFrequencyReductionPlan((1, 2))
        )
        reverse = @inferred KC.reduce_causal_frequency_expression(
            expression, KC.CausalFrequencyReductionPlan((2, 1))
        )
        @test KC.causal_frequency_reduction_kind(forward) === KC.CausalFrequencyIntegrated
        @test KC.causal_frequency_reduction_kind(reverse) === KC.CausalFrequencyIntegrated
        @test KC.causal_frequency_reduction_expression(forward) ==
            KC.causal_frequency_reduction_expression(reverse)

        terms = KC.causal_frequency_terms(KC.causal_frequency_reduction_expression(forward))
        @test length(terms) == 1
        expected = KC.CausalFrequencyTerm(
            -1 // 1,
            [causal_expression_denominator((0, 0), Eₐ + Eᵦ - Eᵧ, 3)],
        )
        @test only(terms) == expected
    end

    @testset "reduction plan preserves typed blocker provenance" begin
        marginal = one_term_expression(
            1 // 1, [causal_expression_denominator((1, 0), -Eₐ, 1)]
        )
        plan = KC.CausalFrequencyReductionPlan((1, 2))
        result = @inferred KC.reduce_causal_frequency_expression(marginal, plan)
        @test KC.causal_frequency_reduction_kind(result) ===
            KC.CausalFrequencyMarginalAtInfinity
        @test KC.causal_frequency_reduction_steps(result) == 0
        @test KC.causal_frequency_reduction_blocked_frequency(result) == 1
    end
end
