using KeldyshContraction, Test
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus

@qfields c::Destroy(Classical) q::Destroy(Quantum)

@testset "is regular" begin
    using KeldyshContraction: regular
    @test !regular((q(Plus), c'))
    @test regular((q(Minus), c'))
    @test !regular((c(Minus), q'))
    @test regular((c(Plus), q'))
    @test regular((c, q'))
end

@testset "regularise vs no regularise in elastic setting" begin
    using KeldyshContraction
    using KeldyshContraction: set_reg_to_zero

    elasctic2boson_reguralisation =
        -0.5 * ((c(Minus)^2 + q(Minus)^2) * c' * q' + c(Plus) * q(Plus) * ((c')^2 + (q')^2))
    elasctic2boson = -0.5 * ((c^2 + q^2) * c' * q' + c * q * ((c')^2 + (q')^2))

    L_reg = InteractionLagrangian(elasctic2boson_reguralisation)
    L_no_reg = InteractionLagrangian(elasctic2boson)

    GF_reg = DressedPropagator(L_reg; simplify=true, _set_reg_to_zero=true)
    GF_no_reg = DressedPropagator(L_no_reg; simplify=true)

    @test isequal(GF_reg.keldysh, GF_no_reg.keldysh)
    @test isequal(GF_reg.retarded, GF_no_reg.retarded)
    @test isequal(GF_reg.advanced, GF_no_reg.advanced)

    GF_reg2 = DressedPropagator(L_reg, 2; simplify=true, _set_reg_to_zero=true)
    GF_no_reg2 = DressedPropagator(L_no_reg, 2; simplify=true)

    @test isequal(GF_reg2.keldysh, GF_no_reg2.keldysh)
    @test isequal(GF_reg2.retarded, GF_no_reg2.retarded)
    @test isequal(GF_reg2.advanced, GF_no_reg2.advanced)
end

@testset "regularise vs no regularise in inelastic setting" begin
    using KeldyshContraction

    loss2boson_unregular =
        im * (
            0.5 * c' * q' * (c^2 + q^2) - 0.5 * c * q * ((c')^2 + (q')^2) +
            c' * q' * (c * q + c * q)
        )

    loss2boson =
        im * (
            0.5 * c' * q' * (c(Minus) * c(Minus) + q(Minus) * q(Minus)) -
            0.5 * c(Plus) * q(Plus) * (c' * c' + q' * q') +
            c' * q' * (c(Plus) * q(Plus) + c(Minus) * q(Minus))
        )

    L_reg = InteractionLagrangian(loss2boson)
    L_no_reg = InteractionLagrangian(loss2boson_unregular)

    GF_reg2 = DressedPropagator(L_reg, 2; simplify=true, _set_reg_to_zero=true)
    GF_no_reg2 = DressedPropagator(L_no_reg, 2; simplify=true, _set_reg_to_zero=true)

    topology_reg = topologies(GF_reg2.keldysh)
    topology_no_reg = topologies(GF_no_reg2.keldysh)

    # The [3] multiplicity doesn't contain equal-time propagators, so the reguralisation shoudn't have any affect.
    @test isequal(Set(topology_reg[[3]]), Set(topology_no_reg[[3]]))
    @test !isequal(Set(topology_reg[[2]]), Set(topology_no_reg[[2]]))
    @test !isequal(Set(topology_reg[[1]]), Set(topology_no_reg[[1]]))
end
