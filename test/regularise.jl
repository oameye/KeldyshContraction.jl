using KeldyshContraction, Test
using KeldyshContraction: Classical, Quantum, Plus, Minus

@qfields c::Destroy(Classical) q::Destroy(Quantum)

@testset "is regular" begin
    using KeldyshContraction: regular
    @test !regular((q(Plus), c'))
    @test regular((q(Minus), c'))
    @test !regular((c(Minus), q'))
    @test regular((c(Plus), q'))
    @test regular((c, q'))
end

@testset "regularise vs no regularise" begin
    using KeldyshContraction

    elasctic2boson_reguralisation =
        -0.5 * ((c(Minus)^2 + q(Minus)^2) * c' * q' + c(Plus) * q(Plus) * ((c')^2 + (q')^2))
    elasctic2boson = -0.5 * ((c^2 + q^2) * c' * q' + c * q * ((c')^2 + (q')^2))

    L_reg = InteractionLagrangian(elasctic2boson_reguralisation)
    L_no_reg = InteractionLagrangian(elasctic2boson)

    GF_reg = DressedPropagator(L_reg)
    GF_no_reg = DressedPropagator(L_no_reg)

    @test isequal(GF_reg.keldysh, GF_no_reg.keldysh)
    @test isequal(GF_reg.retarded, GF_no_reg.retarded)
    @test isequal(GF_reg.advanced, GF_no_reg.advanced)

    GF_reg2 = DressedPropagator(L_reg, 2)
    GF_no_reg2 = DressedPropagator(L_no_reg, 2)

    @test isequal(GF_reg2.keldysh, GF_no_reg2.keldysh)
    @test isequal(GF_reg2.retarded, GF_no_reg2.retarded)
    @test isequal(GF_reg2.advanced, GF_no_reg2.advanced)
end
