function benchmark_two_body_scattering!(SUITE)
    @qfields c::Destroy(Classical) q::Destroy(Quantum)
    elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
    L_int = InteractionLagrangian(elasctic2boson)

    GF = DressedPropagator(L_int)
    # Σ = SelfEnergy(GF)

    SUITE["Two body scattering"]["Green's function"] = @benchmarkable DressedPropagator(
        $L_int;
    ) seconds = 10
    # SUITE["Two body loss"]["Self-energy"] = @benchmarkable SelfEnergy($GF) seconds = 10

    order = 2
    SUITE["Two body scattering"]["Green's function second order"] = @benchmarkable DressedPropagator(
        $L_int; order=($order)
    ) seconds = 50
    return nothing
end
