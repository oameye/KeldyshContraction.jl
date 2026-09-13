import KeldyshContraction as KC

@qfields benchmark_causal_ϕ::Boson

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

    suite["Frequency reduction"]["causal expression residue"] = @benchmarkable KC.integrate_causal_frequency_expression(
        $expression, 1
    ) seconds = 10
    suite["Frequency reduction"]["causal two-frequency plan"] = @benchmarkable KC.reduce_causal_frequency_expression(
        $expression, $plan
    ) seconds = 10
    return suite
end
