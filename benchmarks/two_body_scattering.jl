function benchmark_two_body_scattering!(SUITE)
    @qfields c::Destroy(Classical) q::Destroy(Quantum)
    elasctic2boson = 0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2)
    L_int = InteractionLagrangian(elasctic2boson)

    GF = wick_contraction(L_int)
    # Σ = SelfEnergy(GF)

    simplify = false
    SUITE["Two body scattering"]["Green's function"] = @benchmarkable wick_contraction(
        $L_int; simplify=($(simplify))
    ) seconds = 10
    # SUITE["Two body loss"]["Self-energy"] = @benchmarkable SelfEnergy($GF) seconds = 10

    order = 2
    SUITE["Two body scattering"]["Green's function second order"] = @benchmarkable wick_contraction(
        $L_int; order=($order), simplify=($(simplify))
    ) seconds = 50
    return nothing
end
