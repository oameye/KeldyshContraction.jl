using KeldyshContraction: KeldyshContraction.Regularisation.Plus as Plus
using KeldyshContraction: KeldyshContraction.Regularisation.Minus as Minus

function benchmark_two_body_loss!(SUITE)
    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)

    loss2boson =
        0.5 * ϕ' * ψ' * (ϕ(Minus) * ϕ(Minus) + ψ(Minus) * ψ(Minus)) -
        0.5 * ϕ(Plus) * ψ(Plus) * (ϕ' * ϕ' + ψ' * ψ') +
        ϕ' * ψ' * (ϕ(Plus) * ψ(Plus) + ϕ(Minus) * ψ(Minus))
    L_int = InteractionLagrangian(loss2boson)

    GF = DressedPropagator(L_int)
    Σ = SelfEnergy(GF)

    simplify = true
    _set_reg_to_zero=true
    SUITE["Two body loss"]["Green's function"] = @benchmarkable DressedPropagator($L_int;simplify,_set_reg_to_zero) seconds =
        10
    SUITE["Two body loss"]["Self-energy"] = @benchmarkable SelfEnergy($GF) seconds = 10

    order = 2
    SUITE["Two body loss"]["Green's function second order"] = @benchmarkable DressedPropagator(
        $L_int, $order;simplify,_set_reg_to_zero
    ) seconds = 50
    return nothing
end
