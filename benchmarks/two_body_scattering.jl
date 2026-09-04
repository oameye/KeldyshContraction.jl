function benchmark_two_body_scattering!(SUITE)
    @qfields c::Destroy(Classical) q::Destroy(Quantum)
    elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
    L_int = InteractionLagrangian(elasctic2boson)

    GF = DressedPropagator(L_int, Val(1), Val(3))
    # Σ = SelfEnergy(GF, Val(1))

    SUITE["Two body scattering"]["Green's function"] = @benchmarkable DressedPropagator(
        $L_int, Val(1), Val(3);
    ) seconds = 10
    # SUITE["Two body loss"]["Self-energy"] = @benchmarkable SelfEnergy($GF, Val(1)) seconds = 10

    order = 2
    SUITE["Two body scattering"]["Green's function second order"] = @benchmarkable DressedPropagator(
        $L_int, Val($order), Val(5)
    ) seconds = 50

    GF2 = DressedPropagator(L_int, Val(2), Val(5))
    wigner_transform(GF2)

    SUITE["Two body scattering"]["Wigner transform"] = @benchmarkable wigner_transform($GF2) seconds =
        10

    Σ2 = SelfEnergy(GF2, Val(2))
    Σk2 = wigner_transform(Σ2)
    SUITE["Two body scattering"]["Collision integral"] = @benchmarkable KeldyshContraction.CollisionIntegral(
        $Σk2
    ) seconds = 10
    return nothing
end
