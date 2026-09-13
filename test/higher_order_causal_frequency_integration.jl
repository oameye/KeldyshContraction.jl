using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields higher_order_causal_ϕ::Boson

function higher_order_causal_energies()
    basis = KC.MomentumBasis(2)
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    return EnergyForm(DispersionAtom(higher_order_causal_ϕ, k)),
    EnergyForm(DispersionAtom(higher_order_causal_ϕ, q))
end

function higher_order_denominator(coefficient, energy, infinitesimal)
    return KC.CausalFrequencyDenominator([coefficient // 1], energy, infinitesimal // 1)
end

function higher_order_expression(coefficient, denominators)
    return KC.CausalFrequencyExpression([KC.CausalFrequencyTerm(coefficient, denominators)])
end

function integrated_higher_order_expression(expression)
    result = @inferred KC.integrate_causal_frequency_expression(expression, 1)
    @test KC.causal_frequency_integration_kind(result) === KC.CausalFrequencyIntegrated
    return KC.causal_frequency_integration_expression(result)
end

@testset "higher-order causal frequency integration" begin
    Eₖ, E_q = higher_order_causal_energies()

    @testset "normalized proportional factors form one exact pole class" begin
        first = higher_order_denominator(2, -2 * E_q, -2)
        proportional = higher_order_denominator(-3, 3 * E_q, 3)
        term = KC.CausalFrequencyTerm(1 // 1, [first, proportional])
        denominators = KC.causal_frequency_denominators(term)

        @test length(denominators) == 2
        @test isequal(denominators[1], denominators[2])
        @test KC.causal_frequency_coefficient(term) == -1 // 6
        @test isempty(
            integrated_higher_order_expression(KC.CausalFrequencyExpression([term]))
        )
    end

    @testset "simple pole path is unchanged" begin
        upper = higher_order_denominator(1, -E_q, -1)
        lower = higher_order_denominator(1, -Eₖ, 1)
        expression = higher_order_expression(3 // 2, [upper, lower])

        integrated = integrated_higher_order_expression(expression)
        transformed = higher_order_denominator(0, E_q - Eₖ, 2)
        expected = KC.CausalFrequencyExpression([
            KC.CausalFrequencyTerm(
                convert(KC.ComplexRationals, (3 // 2) * im), [transformed]
            ),
        ])
        @test integrated == expected
    end

    @testset "pure same-side double pole integrates to zero" begin
        upper = higher_order_denominator(1, -E_q, -1)
        expression = higher_order_expression(1 // 1, [upper, upper])
        @test isempty(integrated_higher_order_expression(expression))
    end

    @testset "double-pole derivative residue" begin
        upper = higher_order_denominator(1, -E_q, -1)
        lower = higher_order_denominator(1, -Eₖ, 1)
        expression = higher_order_expression(1 // 1, [upper, upper, lower])

        integrated = integrated_higher_order_expression(expression)
        transformed = higher_order_denominator(0, E_q - Eₖ, 2)
        expected = KC.CausalFrequencyExpression([
            KC.CausalFrequencyTerm(
                convert(KC.ComplexRationals, -im), [transformed, transformed]
            ),
        ])
        @test integrated == expected
    end

    @testset "triple-pole derivative residue includes inverse factorial" begin
        upper = higher_order_denominator(1, -E_q, -1)
        lower = higher_order_denominator(1, -Eₖ, 1)
        expression = higher_order_expression(1 // 1, [upper, upper, upper, lower])

        integrated = integrated_higher_order_expression(expression)
        transformed = higher_order_denominator(0, E_q - Eₖ, 2)
        expected = KC.CausalFrequencyExpression([
            KC.CausalFrequencyTerm(
                convert(KC.ComplexRationals, im), [transformed, transformed, transformed]
            ),
        ])
        @test integrated == expected
    end

    @testset "non-unit proportional coefficients preserve exact Jacobian sign" begin
        first = higher_order_denominator(2, -2 * E_q, -2)
        second = higher_order_denominator(-3, 3 * E_q, 3)
        lower = higher_order_denominator(5, -5 * Eₖ, 5)
        expression = higher_order_expression(1 // 1, [first, second, lower])

        integrated = integrated_higher_order_expression(expression)
        transformed = higher_order_denominator(0, 5 * (E_q - Eₖ), 10)
        expected = KC.CausalFrequencyExpression([
            KC.CausalFrequencyTerm(
                convert(KC.ComplexRationals, (5 // 6) * im), [transformed, transformed]
            ),
        ])
        @test integrated == expected
    end

    @testset "multiple upper pole classes contribute once per class" begin
        upper_double = higher_order_denominator(1, -E_q, -1)
        upper_simple = higher_order_denominator(1, -2 * E_q, -1)
        lower = higher_order_denominator(1, -Eₖ, 1)
        expression = higher_order_expression(
            1 // 1, [upper_double, upper_double, upper_simple, lower]
        )

        integrated = integrated_higher_order_expression(expression)
        simple_at_double = higher_order_denominator(0, -E_q, 0)
        lower_at_double = higher_order_denominator(0, E_q - Eₖ, 2)
        double_at_simple = higher_order_denominator(0, E_q, 0)
        lower_at_simple = higher_order_denominator(0, 2 * E_q - Eₖ, 2)
        minus_i = convert(KC.ComplexRationals, -im)
        plus_i = convert(KC.ComplexRationals, im)
        expected = KC.CausalFrequencyExpression([
            KC.CausalFrequencyTerm(
                minus_i, [simple_at_double, simple_at_double, lower_at_double]
            ),
            KC.CausalFrequencyTerm(
                minus_i, [simple_at_double, lower_at_double, lower_at_double]
            ),
            KC.CausalFrequencyTerm(
                plus_i, [double_at_simple, double_at_simple, lower_at_simple]
            ),
        ])
        @test integrated == expected
    end

    @testset "repeated differentiation combines identical monomials" begin
        upper = higher_order_denominator(1, -E_q, -1)
        lower = higher_order_denominator(1, -Eₖ, 1)
        expression = higher_order_expression(1 // 1, [upper, upper, lower, lower])

        integrated = integrated_higher_order_expression(expression)
        transformed = higher_order_denominator(0, E_q - Eₖ, 2)
        expected = KC.CausalFrequencyExpression([
            KC.CausalFrequencyTerm(
                convert(KC.ComplexRationals, -2im), [transformed, transformed, transformed]
            ),
        ])
        @test integrated == expected
    end

    @testset "denominator permutation is canonical" begin
        upper = higher_order_denominator(1, -E_q, -1)
        lower = higher_order_denominator(1, -Eₖ, 1)
        first = higher_order_expression(2 // 3, [upper, upper, lower])
        second = higher_order_expression(2 // 3, [lower, upper, upper])

        @test integrated_higher_order_expression(first) ==
            integrated_higher_order_expression(second)
    end

    @testset "opposite prescriptions remain a pinch" begin
        lower = higher_order_denominator(1, -E_q, 1)
        upper = higher_order_denominator(1, -E_q, -1)
        expression = higher_order_expression(1 // 1, [lower, upper])

        canonical_term = only(KC.causal_frequency_terms(expression))
        denominators = KC.causal_frequency_denominators(canonical_term)
        @test !isequal(denominators[1], denominators[2])

        result = @inferred KC.integrate_causal_frequency_expression(expression, 1)
        @test KC.causal_frequency_integration_kind(result) === KC.CausalFrequencyPinch
        @test KC.causal_frequency_integration_expression(result) ==
            KC._complex_causal_frequency_expression(expression)
    end
end
