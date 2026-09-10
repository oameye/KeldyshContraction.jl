using KeldyshContraction

@qfields benchmark_collision_ϕ::Boson

function benchmark_off_shell_collision!(SUITE)
    c = benchmark_collision_ϕ[Classical]
    q = benchmark_collision_ϕ[Quantum]
    interaction =
        -(1 // 2) * ((c^2 + q^2) * bar(c) * bar(q) + c * q * (bar(c)^2 + bar(q)^2))
    L = InteractionLagrangian(interaction, :g)
    G = DressedPropagator(L, Val(2), Val(5); simplify=false)
    ΣW = wigner_transform(SelfEnergy(fourier_transform(G)); gradient_order=Val(0))
    KΣ = kinetic_expression(ΣW)

    SUITE["Off-shell collision"]["Kadanoff-Baym identity"] = @benchmarkable off_shell_collision_expression(
        $KΣ
    ) seconds = 10
    return nothing
end
