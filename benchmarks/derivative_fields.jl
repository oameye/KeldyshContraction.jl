using KeldyshContraction

function benchmark_derivative_fields!(SUITE)
    @qfields ψderivbench::Fermion
    ψ₁, ψ₂ = ψderivbench[One], ψderivbench[Two]
    ∂xψ₂ = partial(ψ₂, :x)
    vertex = ψ₁ * ∂xψ₂ * bar(ψ₁) * bar(∂xψ₂)
    L = InteractionLagrangian(vertex, :γ)
    G = DressedPropagator(L, Val(1), Val(3); simplify=false)

    SUITE["Derivative fields"]["decorate"] = @benchmarkable partial($ψ₂, :x) seconds = 10
    SUITE["Derivative fields"]["Grassmann product"] = @benchmarkable $ψ₁ * $∂xψ₂ seconds =
        10
    SUITE["Derivative fields"]["propagator order 1"] = @benchmarkable DressedPropagator(
        $L, Val(1), Val(3); simplify=false
    ) seconds = 10
    SUITE["Derivative fields"]["self-energy order 1"] = @benchmarkable SelfEnergy($G) seconds =
        10
    return nothing
end
