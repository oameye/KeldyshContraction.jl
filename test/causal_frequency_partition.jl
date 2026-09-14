using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields causal_partition_ϕ::Boson

function partition_denominator(coefficient, energy, infinitesimal)
    return KC.CausalFrequencyDenominator([coefficient // 1], energy, infinitesimal // 1)
end

@testset "blocked causal stage preserves independent regular subblock" begin
    basis = KC.MomentumBasis(2)
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    Eₖ = EnergyForm(DispersionAtom(causal_partition_ϕ, k))
    E_q = EnergyForm(DispersionAtom(causal_partition_ϕ, q))

    upper_k = partition_denominator(1, -Eₖ, -1)
    lower_k = partition_denominator(1, -Eₖ, 1)
    lower_q = partition_denominator(1, -E_q, 1)

    regular_term = KC.CausalFrequencyTerm(1 // 1, [upper_k, upper_k, lower_q])
    pinch_term = KC.CausalFrequencyTerm(2 // 1, [upper_k, lower_k])
    expression = KC.CausalFrequencyExpression([regular_term, pinch_term])

    complete = @inferred KC.integrate_causal_frequency_expression(expression, 1)
    @test KC.causal_frequency_integration_kind(complete) === KC.CausalFrequencyPinch

    regular, residual = @inferred KC._split_contour_safe_frequency_block(expression, 1)
    @test !isempty(regular)
    @test !isempty(residual)
    @test KC._merge_trotter_causal_expression(regular, residual) ==
        KC._complex_causal_frequency_expression(expression)

    regular_result = @inferred KC.integrate_causal_frequency_expression(regular, 1)
    @test KC.causal_frequency_integration_kind(regular_result) ===
        KC.CausalFrequencyIntegrated
    @test !isempty(KC.causal_frequency_integration_expression(regular_result))

    residual_result = @inferred KC.integrate_causal_frequency_expression(residual, 1)
    @test KC.causal_frequency_integration_kind(residual_result) === KC.CausalFrequencyPinch
    @test KC.causal_frequency_integration_expression(residual_result) == residual
end
