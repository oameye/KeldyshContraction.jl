using KeldyshContraction, Test
import KeldyshContraction as KC

@testset "native perturbation parameters" begin
    @syms g λ
    gₘ = parameter_monomial(g)
    λₘ = parameter_monomial(λ)

    @test gₘ isa ParameterMonomial
    @test @inferred(gₘ * λₘ) == @inferred(λₘ * gₘ)
    @test hash(gₘ * λₘ) == hash(λₘ * gₘ)
    @test parameter_monomial(g * λ) == gₘ * λₘ
    @test parameter_monomial(λ * g) == gₘ * λₘ
    @test parameter_monomial(g^2) == @inferred(gₘ^2)
    @test @inferred(gₘ^0) == one(ParameterMonomial)
    @test eltype(gₘ.powers) === KC.ParameterPower
    @test isconcretetype(typeof(gₘ))
    @test isconcretetype(eltype(gₘ.powers))
end
