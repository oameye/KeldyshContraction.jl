using KeldyshContraction

@qfields benchmark_wigner_ψ::Fermion
const benchmark_wigner_ψ₁ = benchmark_wigner_ψ[One]
const benchmark_wigner_ψ₂ = benchmark_wigner_ψ[Two]

function benchmark_wigner!(SUITE)
    ∂xψ₂ = partial(benchmark_wigner_ψ₂, :x)
    vertex = benchmark_wigner_ψ₁ * ∂xψ₂ * bar(benchmark_wigner_ψ₁) * bar(∂xψ₂)
    L = InteractionLagrangian(vertex, :d)
    G = DressedPropagator(L, Val(1), Val(3); target=benchmark_wigner_ψ, simplify=false)
    Gk = fourier_transform(G)
    Σk = SelfEnergy(Gk)

    SUITE["Wigner representation"]["p-wave propagator order 0"] = @benchmarkable wigner_transform(
        $Gk; gradient_order=Val(0)
    ) seconds = 10
    SUITE["Wigner representation"]["p-wave self-energy order 0"] = @benchmarkable wigner_transform(
        $Σk; gradient_order=Val(0)
    ) seconds = 10
    return nothing
end
