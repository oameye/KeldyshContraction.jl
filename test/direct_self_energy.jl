using KeldyshContraction, Test
import KeldyshContraction as KC

function test_direct_self_energy_equal(direct, reference)
    @test direct.keldysh == reference.keldysh
    @test direct.retarded == reference.retarded
    @test direct.advanced == reference.advanced
    @test KC.parameters(direct) == KC.parameters(reference)
    @test KC.target_family(direct) == KC.target_family(reference)
    @test KC.order(direct) == KC.order(reference)
end

@testset "Direct self-energy generation" begin
    @qfields direct_Σ_ϕ::Boson
    c, q = direct_Σ_ϕ[Classical], direct_Σ_ϕ[Quantum]
    elastic = -(0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2))
    L = InteractionLagrangian(elastic, :g)

    for (order, edges) in ((Val(2), Val(5)), (Val(3), Val(7)))
        reference = SelfEnergy(DressedPropagator(L, order, edges))
        direct = SelfEnergy(L, order, edges)
        test_direct_self_energy_equal(direct, reference)
    end

    loss =
        0.5 *
        bar(c) *
        bar(q) *
        (
            c(KC.Regularisation.Minus) * c(KC.Regularisation.Minus) +
            q(KC.Regularisation.Minus) * q(KC.Regularisation.Minus)
        ) -
        0.5 *
        c(KC.Regularisation.Plus) *
        q(KC.Regularisation.Plus) *
        (bar(c) * bar(c) + bar(q) * bar(q)) +
        bar(c) *
        bar(q) *
        (
            c(KC.Regularisation.Plus) * q(KC.Regularisation.Plus) +
            c(KC.Regularisation.Minus) * q(KC.Regularisation.Minus)
        )
    L_loss = InteractionLagrangian(loss, :γ)
    reference_loss = SelfEnergy(
        DressedPropagator(
            L_loss, Val(2), Val(5); preserve_regularisation=true
        )
    )
    direct_loss = SelfEnergy(
        L_loss, Val(2), Val(5); preserve_regularisation=true
    )
    test_direct_self_energy_equal(direct_loss, reference_loss)
end
