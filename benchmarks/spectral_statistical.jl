using KeldyshContraction

@qfields benchmark_kinetic_ψ::Fermion
const benchmark_kinetic_ψ₁ = benchmark_kinetic_ψ[One]
const benchmark_kinetic_ψ₂ = benchmark_kinetic_ψ[Two]

function benchmark_spectral_statistical!(SUITE)
    ∂xψ₂ = partial(benchmark_kinetic_ψ₂, :x)
    vertex = benchmark_kinetic_ψ₁ * ∂xψ₂ * bar(benchmark_kinetic_ψ₁) * bar(∂xψ₂)
    L = InteractionLagrangian(vertex, :d)
    G = DressedPropagator(L, Val(1), Val(3); target=benchmark_kinetic_ψ, simplify=false)
    ΣW = wigner_transform(SelfEnergy(fourier_transform(G)); gradient_order=Val(0))
    KΣ = kinetic_expression(ΣW)

    SUITE["Spectral statistical IR"]["p-wave self-energy lowering"] = @benchmarkable kinetic_expression(
        $ΣW
    ) seconds = 10
    SUITE["Spectral statistical IR"]["retarded minus advanced"] = @benchmarkable retarded_minus_advanced(
        $KΣ
    ) seconds = 10
    return nothing
end
