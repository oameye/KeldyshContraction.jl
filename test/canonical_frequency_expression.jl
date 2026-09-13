using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields canonical_frequency_φ::Boson

function canonical_frequency_energy(family, momentum)
    return KC.EnergyForm(KC.DispersionAtom(family, momentum))
end

@testset "Canonical causal frequency expression" begin
    basis = KC.MomentumBasis(3)
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    E_k = canonical_frequency_energy(canonical_frequency_φ, k)
    E_q = canonical_frequency_energy(canonical_frequency_φ, q)
    E_r = canonical_frequency_energy(canonical_frequency_φ, r)

    denominator = KC.CausalFrequencyDenominator([1 // 1], -E_k, 1 // 1)
    scaled = KC.CausalFrequencyDenominator([-2 // 1], 2 * E_k, -2 // 1)
    canonical, factor = KC.canonical_causal_frequency_denominator(scaled)
    @test factor == -1 // 2
    @test canonical == denominator

    original_term = KC.CausalFrequencyTerm(1 // 1, [denominator])
    scaled_term = KC.CausalFrequencyTerm(2 // 1, [scaled])
    @test KC.CausalFrequencyExpression([original_term, scaled_term]) |> isempty

    # Each term separately has a marginal 1/ω tail, but the complete expression does not:
    #
    #   1/(ω-E_k+i0) - 1/(ω-E_q+i0) = O(ω^-2).
    #
    # Reducibility therefore cannot be decided term by term.
    d_k = KC.CausalFrequencyDenominator([1 // 1], -E_k, 1 // 1)
    d_q = KC.CausalFrequencyDenominator([1 // 1], -E_q, 1 // 1)
    marginal = KC.CausalFrequencyExpression([
        KC.CausalFrequencyTerm(1 // 1, [d_k]), KC.CausalFrequencyTerm(-1 // 1, [d_q])
    ])
    @test KC.frequency_decay_lower_bound(marginal, 1) == 2
    @test KC.causal_contour_safe_at_infinity(marginal, 1)

    single = KC.CausalFrequencyExpression([KC.CausalFrequencyTerm(1 // 1, [d_k])])
    @test KC.frequency_decay_lower_bound(single, 1) == 1
    @test !KC.causal_contour_safe_at_infinity(single, 1)

    # The same cancellation must be visible with a nontrivial spectator rational factor.
    d_k_2d = KC.CausalFrequencyDenominator([1 // 1, 0 // 1], -E_k, 1 // 1)
    d_q_2d = KC.CausalFrequencyDenominator([1 // 1, 0 // 1], -E_q, 1 // 1)
    spectator = KC.CausalFrequencyDenominator([0 // 1, 1 // 1], -E_r, -1 // 1)
    with_spectator = KC.CausalFrequencyExpression([
        KC.CausalFrequencyTerm(1 // 1, [d_k_2d, spectator]),
        KC.CausalFrequencyTerm(-1 // 1, [d_q_2d, spectator]),
    ])
    order, leading = KC.leading_frequency_asymptotic(with_spectator, 1)
    @test order == 1
    @test isempty(leading)
    @test KC.causal_contour_safe_at_infinity(with_spectator, 1)

    @test @inferred(KC.CausalFrequencyTerm(1 // 1, [denominator])) == original_term
end
