using KeldyshContraction, Test
using KeldyshContraction: In, Out
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus
using KeldyshContraction: Bulk, position

@qfields ϕᶜ::Destroy(Classical) ϕᴾ::Destroy(Quantum)

L_int =
    im*(
        0.5 * ϕᶜ' * ϕᴾ' * (ϕᶜ(Minus) * ϕᶜ(Minus) + ϕᴾ(Minus) * ϕᴾ(Minus)) -
        0.5 * ϕᶜ(Plus) * ϕᴾ(Plus) * (ϕᶜ' * ϕᶜ' + ϕᴾ' * ϕᴾ') +
        ϕᶜ' * ϕᴾ' * (ϕᶜ(Plus) * ϕᴾ(Plus) + ϕᶜ(Minus) * ϕᴾ(Minus))
    )
L = InteractionLagrangian(L_int)

@testset "set_position" begin
    @test position(L) == Bulk(1)
    @test position(L(2)) == Bulk(2)
end

@testset "ladder sorted" begin
    L1 = L
    L2 = L(2)
    expr = L1.lagrangian * L2.lagrangian
    for arg in expr.arguments
        sorted = sort(arg.args_nc; by=KeldyshContraction.ladder)
        @test isequal(arg.args_nc, sorted)
    end
end

@testset "zero loop filter" begin
    using KeldyshContraction: has_zero_loop, Contraction
    # Gᴿ(1,2) Gᴿ(2,1) = 0
    vs = Contraction[(ϕᶜ(Bulk(1)), ϕᴾ'(Bulk(2))), (ϕᶜ(Bulk(2)), ϕᴾ'(Bulk(1)))]
    @test has_zero_loop(vs)

    # Gᴬ(1,2) Gᴬ(2,1) = 0
    vs = Contraction[(ϕᴾ(Bulk(1)), ϕᶜ'(Bulk(2))), (ϕᴾ(Bulk(2)), ϕᶜ'(Bulk(1)))]
    @test has_zero_loop(vs)

    # Gᴿ(1,2) Gᴬ(1,2) = 0
    vs = Contraction[(ϕᶜ(Bulk(1)), ϕᴾ'(Bulk(2))), (ϕᴾ(Bulk(1)), ϕᶜ'(Bulk(2)))]
    @test has_zero_loop(vs)

    vs = Contraction[(ϕᶜ(Bulk(1)), ϕᴾ'(Bulk(2))), (ϕᶜ(Bulk(1)), ϕᴾ'(Bulk(2)))]
    @test !has_zero_loop(vs)
end

@testset "canonicalize" begin
    using KeldyshContraction: canonicalize, Contraction, is_canonical, Out, In
    vs = Contraction[(ϕᶜ(Out()), ϕᴾ'(Bulk(2))), (ϕᶜ(Bulk(2)), ϕᴾ'(Bulk(1)))]
    @test !is_canonical(vs)

    vs′ = Contraction[(ϕᶜ(Out()), ϕᴾ'(Bulk(1))), (ϕᶜ(Bulk(1)), ϕᴾ'(Bulk(2)))]
    @test is_canonical(vs′)

    @test canonicalize(vs) == vs′
end
