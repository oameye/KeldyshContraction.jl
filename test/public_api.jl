using KeldyshContraction, Test
import KeldyshContraction as KC

@testset "declared public API" begin
    exported = Set(names(KC; all=false, imported=false))
    @test :InteractionLagrangian in exported
    @test :DressedPropagator in exported
    @test :SelfEnergy in exported
    @test :fourier_transform in exported
    @test :wigner_transform in exported
    @test :kinetic_expression in exported
    @test :collision_kernel in exported

    if VERSION >= v"1.11"
        @test Base.ispublic(KC, :QTerm)
        @test Base.ispublic(KC, :Regularisation)
        @test Base.ispublic(KC, :MomentumBasis)
        @test Base.ispublic(KC, :matrix)
    end

    @test !isdefined(KC, :Destroy)
    @test !isdefined(KC, :Create)
end
