import KeldyshContraction as KC

@qfields benchmark_causal_ϕ::Boson benchmark_trotter_ϕ::Boson

function benchmark_trotter_collision()
    c = benchmark_trotter_ϕ[Classical]
    q = benchmark_trotter_ϕ[Quantum]
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus
    loss =
        (1 // 2) *
        im *
        (
            bar(c) * bar(q) * (c(minus)^2 + q(minus)^2) -
            c(plus) * q(plus) * (bar(c)^2 + bar(q)^2) +
            2 * bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
        )

    L = InteractionLagrangian(loss, :γ)
    G = DressedPropagator(L, Val(1), Val(3); simplify=true, _set_reg_to_zero=false)
    GF = fourier_transform(G)
    ΣF = SelfEnergy(GF)
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    KΣ = kinetic_expression(ΣW)
    off_shell = off_shell_collision_expression(KΣ)
    spectral = spectral_dispersive_collision(off_shell)
    return KC.canonical_frequency_collision(spectral)
end

function benchmark_causal_frequency_integration!(suite)
    basis = KC.MomentumBasis(3)
    momenta = ntuple(i -> KC.basis_momentum(basis, i), 3)
    energies = map(
        momentum -> KC.EnergyForm(KC.DispersionAtom(benchmark_causal_ϕ, momentum)), momenta
    )
    Eₐ, Eᵦ, Eᵧ = energies

    upper_1 = KC.CausalFrequencyDenominator([1 // 1, 0 // 1], -Eₐ, -1 // 1)
    upper_2 = KC.CausalFrequencyDenominator([0 // 1, 1 // 1], -Eᵦ, -1 // 1)
    coupled_lower = KC.CausalFrequencyDenominator([1 // 1, 1 // 1], -Eᵧ, 1 // 1)
    expression = KC.CausalFrequencyExpression([
        KC.CausalFrequencyTerm(1 // 1, [upper_1, upper_2, coupled_lower])
    ])
    plan = KC.CausalFrequencyReductionPlan((1, 2))

    zero_energy = KC.EnergyForm{Boson}(3)
    support_denominators = KC.CausalFrequencyDenominator{Boson}[
        KC.CausalFrequencyDenominator([1 // 1, 0 // 1], zero_energy, -1 // 1),
        KC.CausalFrequencyDenominator([0 // 1, 1 // 1], zero_energy, -1 // 1),
        KC.CausalFrequencyDenominator([1 // 1, 1 // 1], zero_energy, 1 // 1),
    ]
    pinch_expression = KC.CausalFrequencyExpression([
        KC.CausalFrequencyTerm(1 // 1, support_denominators)
    ])

    trotter_collision = benchmark_trotter_collision()
    shifted_term = first(KC.shifted_frequency_terms(trotter_collision))
    trotter_state = KC.TrotterFrequencyState(shifted_term)
    trotter_witness = KC.trotter_frequency_witness(trotter_state)

    suite["Frequency reduction"]["causal expression residue"] = @benchmarkable KC.integrate_causal_frequency_expression(
        $expression, 1
    ) seconds = 10
    suite["Frequency reduction"]["causal two-frequency plan"] = @benchmarkable KC.reduce_causal_frequency_expression(
        $expression, $plan
    ) seconds = 10
    suite["Frequency reduction"]["affine singular support"] = @benchmarkable KC.affine_singular_support(
        $support_denominators
    ) seconds = 10
    suite["Frequency reduction"]["support-aware pinch plan"] = @benchmarkable KC.reduce_causal_frequency_with_support(
        $pinch_expression, $plan
    ) seconds = 10
    suite["Frequency reduction"]["isolated Trotter elimination"] = @benchmarkable KC.eliminate_trotter_frequency(
        $trotter_state, $trotter_witness
    ) seconds = 10
    suite["Frequency reduction"]["group shifted Trotter collision"] = @benchmarkable KC.reduce_shifted_trotter_frequencies(
        $trotter_collision
    ) seconds = 10
    return suite
end
