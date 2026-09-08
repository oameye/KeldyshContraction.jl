using KeldyshContraction, Test
import KeldyshContraction as KC
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus

@qfields probe_ϕ::Boson
const probe_c = probe_ϕ[Classical]
const probe_q = probe_ϕ[Quantum]
@syms probe_Γ probe_g

function print_collision(label, Σ)
    println("=== ", label, " SELF ENERGY ===")
    println("K=", repr(Σ.keldysh))
    println("R=", repr(Σ.retarded))
    println("A=", repr(Σ.advanced))
    try
        ci = KC.CollisionIntegral(Σ)
        println("=== ", label, " COLLISION ===")
        println(repr(ci))
        for topology in sort!(collect(keys(ci.terms)); by=collect)
            println(label, " topology ", collect(topology), " => ", repr(ci.terms[topology]))
        end
    catch err
        println("=== ", label, " COLLISION ERROR ===")
        showerror(stdout, err, catch_backtrace())
        println()
    end
    return nothing
end

@testset "bosonic collision analytic probe" begin
    elastic = -(
        1 // 2 * (probe_c^2 + probe_q^2) * bar(probe_c) * bar(probe_q) +
        1 // 2 * probe_c * probe_q * (bar(probe_c)^2 + bar(probe_q)^2)
    )
    L_elastic = InteractionLagrangian(elastic)
    GF_elastic = DressedPropagator(L_elastic, Val(2), Val(5); simplify=true)
    Σ_elastic = SelfEnergy(GF_elastic)
    print_collision("ELASTIC_G2", Σ_elastic)

    loss = im * (
        1 // 2 * bar(probe_c) * bar(probe_q) *
        (probe_c(Minus) * probe_c(Minus) + probe_q(Minus) * probe_q(Minus)) -
        1 // 2 * probe_c(Plus) * probe_q(Plus) *
        (bar(probe_c) * bar(probe_c) + bar(probe_q) * bar(probe_q)) +
        bar(probe_c) * bar(probe_q) *
        (probe_c(Plus) * probe_q(Plus) + probe_c(Minus) * probe_q(Minus))
    )
    L_loss = InteractionLagrangian(loss)

    GF_loss1 = DressedPropagator(
        L_loss, Val(1), Val(3); simplify=true, _set_reg_to_zero=true
    )
    Σ_loss1 = SelfEnergy(GF_loss1)
    print_collision("LOSS_GAMMA1", Σ_loss1)

    GF_loss2 = DressedPropagator(
        L_loss, Val(2), Val(5); simplify=true, _set_reg_to_zero=true
    )
    Σ_loss2 = SelfEnergy(GF_loss2)
    print_collision("LOSS_GAMMA2", Σ_loss2)

    L_elastic_p = InteractionLagrangian(elastic, probe_g)
    L_loss_p = InteractionLagrangian(loss, probe_Γ)
    GF_sum2 = DressedPropagator(
        L_elastic_p + L_loss_p,
        Val(2),
        Val(5);
        simplify=true,
        _set_reg_to_zero=true,
    )
    Σ_sum2 = SelfEnergy(GF_sum2)
    println("=== MIXED PARAMETERS ===")
    println(repr(parameters(Σ_sum2)))
    print_collision("MIXED_G_GAMMA", Σ_sum2[probe_g * probe_Γ])

    @test true
end
