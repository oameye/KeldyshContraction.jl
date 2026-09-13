using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields affine_support_ϕ::Boson

zero_energy(n=3) = KC.EnergyForm{Boson}(n)
function constraint(coefficients)
    return KC.AffineFrequencyConstraint(coefficients, zero_energy(length(coefficients)))
end

@testset "affine singular support" begin
    @testset "rank-one repeated support" begin
        support = @inferred KC.affine_singular_support([
            constraint([1 // 1, 0 // 1]), constraint([-2 // 1, 0 // 1])
        ])
        @test KC.has_affine_singular_support(support)
        @test KC.affine_support_rank(support) == 1
        @test KC.affine_support_count(support) == 2
        @test support.dependency_rows == [[1 // 1]]
    end

    @testset "row order and scale invariance" begin
        a = KC.affine_singular_support([
            constraint([1 // 1, 0 // 1]),
            constraint([0 // 1, 1 // 1]),
            constraint([1 // 1, 1 // 1]),
        ])
        b = KC.affine_singular_support([
            constraint([-4 // 1, -4 // 1]),
            constraint([6 // 1, 0 // 1]),
            constraint([0 // 1, -3 // 1]),
        ])
        @test a == b
        @test KC.affine_support_rank(a) == 2
        @test KC.affine_support_count(a) == 3
        @test a.dependency_rows == [[1 // 1, 1 // 1]]
    end

    @testset "independent spectators are coloops" begin
        support = @inferred KC.affine_singular_support([
            constraint([1 // 1, 0 // 1, 0 // 1]),
            constraint([0 // 1, 1 // 1, 0 // 1]),
            constraint([1 // 1, 1 // 1, 0 // 1]),
            constraint([0 // 1, 0 // 1, 1 // 1]),
        ])
        @test KC.affine_support_rank(support) == 2
        @test KC.affine_support_count(support) == 3
        @test support.dependency_rows == [[1 // 1, 1 // 1]]
        @test all(support.independent_constraints) do row
            iszero(row.loop_coefficients[3])
        end
    end

    @testset "energy is part of affine support" begin
        p = KC.LinearMomentum([1, 0, 0])
        q = KC.LinearMomentum([0, 1, 0])
        εp = KC.EnergyForm(KC.DispersionAtom(affine_support_ϕ, p))
        εq = KC.EnergyForm(KC.DispersionAtom(affine_support_ϕ, q))
        support = @inferred KC.affine_singular_support([
            KC.AffineFrequencyConstraint([1 // 1], εp),
            KC.AffineFrequencyConstraint([1 // 1], εq),
        ])
        @test !KC.has_affine_singular_support(support)
        @test KC.affine_support_rank(support) == 0
        @test KC.affine_support_count(support) == 0
    end

    @testset "global contour safety remains expression-level" begin
        p = KC.LinearMomentum([1])
        εp = KC.EnergyForm(KC.DispersionAtom(affine_support_ϕ, p))
        d1 = KC.CausalFrequencyDenominator([1 // 1], zero_energy(1), -1 // 1)
        d2 = KC.CausalFrequencyDenominator([1 // 1], εp, -1 // 1)
        expression = KC.CausalFrequencyExpression([
            KC.CausalFrequencyTerm(1 // 1, [d1]), KC.CausalFrequencyTerm(-1 // 1, [d2])
        ])
        @test KC.frequency_decay_lower_bound(expression, 1) >= 2
        plain = KC.reduce_causal_frequency_expression(
            expression, KC.CausalFrequencyReductionPlan((1,))
        )
        supported = @inferred KC.reduce_causal_frequency_with_support(
            expression, KC.CausalFrequencyReductionPlan((1,))
        )
        @test KC.causal_frequency_support_kind(supported) === KC.CausalFrequencyIntegrated
        @test KC.causal_frequency_support_expression(supported) ==
            KC.causal_frequency_reduction_expression(plain)
        @test isempty(KC.causal_frequency_support_expression(supported))
        @test isempty(KC.causal_frequency_singular_supports(supported))
    end

    @testset "same-side affine dependency is not automatically singular" begin
        d = KC.CausalFrequencyDenominator([1 // 1], zero_energy(1), -1 // 1)
        term = KC.CausalFrequencyTerm(1 // 1, [d, d])
        support = KC.affine_singular_support(term)
        @test KC.has_affine_singular_support(support)
        @test KC.affine_support_rank(support) == 1
        @test KC.affine_support_count(support) == 2

        result = @inferred KC.reduce_causal_frequency_with_support(
            KC.CausalFrequencyExpression([term]), KC.CausalFrequencyReductionPlan((1,))
        )
        @test KC.causal_frequency_support_kind(result) === KC.CausalFrequencyIntegrated
        @test isempty(KC.causal_frequency_support_expression(result))
        @test isempty(KC.causal_frequency_singular_supports(result))
    end

    @testset "support survives pivot-order-dependent pinch" begin
        E = zero_energy(2)
        d1 = KC.CausalFrequencyDenominator([1 // 1, 0 // 1], E, -1 // 1)
        d2 = KC.CausalFrequencyDenominator([0 // 1, 1 // 1], E, -1 // 1)
        d3 = KC.CausalFrequencyDenominator([1 // 1, 1 // 1], E, 1 // 1)
        term = KC.CausalFrequencyTerm(1 // 1, [d1, d2, d3])
        expression = KC.CausalFrequencyExpression([term])

        forward_plan = KC.CausalFrequencyReductionPlan((1, 2))
        reverse_plan = KC.CausalFrequencyReductionPlan((2, 1))
        forward = KC.reduce_causal_frequency_expression(expression, forward_plan)
        reverse = KC.reduce_causal_frequency_expression(expression, reverse_plan)
        @test KC.causal_frequency_reduction_kind(forward) === KC.CausalFrequencyPinch
        @test KC.causal_frequency_reduction_kind(reverse) === KC.CausalFrequencyPinch
        @test KC.causal_frequency_reduction_blocked_frequency(forward) == 2
        @test KC.causal_frequency_reduction_blocked_frequency(reverse) == 1

        support = @inferred KC.affine_singular_support(term)
        @test KC.affine_support_rank(support) == 2
        @test KC.affine_support_count(support) == 3
        @test support.dependency_rows == [[1 // 1, 1 // 1]]

        supported_forward = @inferred KC.reduce_causal_frequency_with_support(
            expression, forward_plan
        )
        supported_reverse = @inferred KC.reduce_causal_frequency_with_support(
            expression, reverse_plan
        )
        @test KC.causal_frequency_support_kind(supported_forward) ===
            KC.CausalFrequencyPinch
        @test KC.causal_frequency_support_kind(supported_reverse) ===
            KC.CausalFrequencyPinch
        @test KC.causal_frequency_support_blocked_frequency(supported_forward) == 2
        @test KC.causal_frequency_support_blocked_frequency(supported_reverse) == 1
        @test KC.causal_frequency_singular_supports(supported_forward) == [support]
        @test KC.causal_frequency_singular_supports(supported_reverse) == [support]
    end
end
