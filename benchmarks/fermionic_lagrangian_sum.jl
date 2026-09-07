using KeldyshContraction

function benchmark_fermionic_lagrangian_sum!(SUITE)
    @qfields ψsumbench::Fermion
    ψ₁, ψ₂ = ψsumbench[One], ψsumbench[Two]
    vertex = ψ₁ * ψ₂ * bar(ψ₁) * bar(ψ₂)
    L_u = InteractionLagrangian(vertex, :u)
    L_v = InteractionLagrangian(vertex, :v)
    Ls = L_u + L_v

    SUITE["Fermionic LagrangianSum"]["single order 1"] = @benchmarkable DressedPropagator(
        $L_u, Val(1), Val(3); simplify=false
    ) seconds = 10
    SUITE["Fermionic LagrangianSum"]["sum order 1"] = @benchmarkable DressedPropagator(
        $Ls, Val(1), Val(3); simplify=false
    ) seconds = 10
    SUITE["Fermionic LagrangianSum"]["sum order 2"] = @benchmarkable DressedPropagator(
        $Ls, Val(2), Val(5); simplify=false
    ) seconds = 20
    return nothing
end
