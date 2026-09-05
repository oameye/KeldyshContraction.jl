using KeldyshContraction
using KeldyshContraction: KeldyshContraction.Regularisation.Plus as Plus
using KeldyshContraction: KeldyshContraction.Regularisation.Minus as Minus

function benchmark_two_body_loss!(SUITE)
    @qfields ϕ::Boson
    c, q = ϕ[Classical], ϕ[Quantum]

    loss2boson =
        0.5 * bar(c) * bar(q) * (c(Minus) * c(Minus) + q(Minus) * q(Minus)) -
        0.5 * c(Plus) * q(Plus) * (bar(c) * bar(c) + bar(q) * bar(q)) +
        bar(c) * bar(q) * (c(Plus) * q(Plus) + c(Minus) * q(Minus))
    L_int = InteractionLagrangian(loss2boson)

    GF = DressedPropagator(L_int, Val(1), Val(3); simplify=true, _set_reg_to_zero=true)
    Σ = SelfEnergy(GF, Val(1))

    SUITE["Two body loss"]["Green's function"] = @benchmarkable DressedPropagator(
        $L_int, Val(1), Val(3); simplify=true, _set_reg_to_zero=true
    ) seconds = 10
    SUITE["Two body loss"]["Self-energy"] = @benchmarkable SelfEnergy($GF, Val(1)) seconds =
        10

    order = 2
    SUITE["Two body loss"]["Green's function second order"] = @benchmarkable DressedPropagator(
        $L_int, Val($order), Val(5); simplify=true, _set_reg_to_zero=true
    ) seconds = 50
    return nothing
end
