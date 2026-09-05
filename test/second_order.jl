using KeldyshContraction, Test
using KeldyshContraction: In, Out
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus
using KeldyshContraction: Bulk, position

@qfields ϕᶜ::Boson(Classical) ϕᴾ::Boson(Quantum)

L_int =
    im * (
        0.5 * bar(ϕᶜ) * bar(ϕᴾ) * (ϕᶜ(Minus) * ϕᶜ(Minus) + ϕᴾ(Minus) * ϕᴾ(Minus)) -
        0.5 * ϕᶜ(Plus) * ϕᴾ(Plus) * (bar(ϕᶜ) * bar(ϕᶜ) + bar(ϕᴾ) * bar(ϕᴾ)) +
        bar(ϕᶜ) * bar(ϕᴾ) * (ϕᶜ(Plus) * ϕᴾ(Plus) + ϕᶜ(Minus) * ϕᴾ(Minus))
    )
L = InteractionLagrangian(L_int)

@testset "set_position" begin
    @test position(L) == Bulk(1)
    @test position(L(2)) == Bulk(2)
end

@testset "orientation sorted" begin
    L1 = L
    L2 = L(2)
    expr = L1.lagrangian * L2.lagrangian
    for arg in expr.arguments
        sorted = sort(arg.args_nc; by=KeldyshContraction.ladder)
        @test all(
            KeldyshContraction.ladder.(arg.args_nc) .== KeldyshContraction.ladder.(sorted)
        )
    end
end

@testset "zero loop filter" begin
    using KeldyshContraction: has_zero_loop, Contraction
    vs = Contraction[(ϕᶜ(Bulk(1)), bar(ϕᴾ)(Bulk(2))), (ϕᶜ(Bulk(2)), bar(ϕᴾ)(Bulk(1)))]
    @test has_zero_loop(vs)

    vs = Contraction[(ϕᴾ(Bulk(1)), bar(ϕᶜ)(Bulk(2))), (ϕᴾ(Bulk(2)), bar(ϕᶜ)(Bulk(1)))]
    @test has_zero_loop(vs)

    vs = Contraction[(ϕᶜ(Bulk(1)), bar(ϕᴾ)(Bulk(2))), (ϕᴾ(Bulk(1)), bar(ϕᶜ)(Bulk(2)))]
    @test has_zero_loop(vs)

    vs = Contraction[(ϕᶜ(Bulk(1)), bar(ϕᴾ)(Bulk(2))), (ϕᶜ(Bulk(1)), bar(ϕᴾ)(Bulk(2)))]
    @test !has_zero_loop(vs)

    vs = Contraction[
        (ϕᶜ(Bulk(1)), bar(ϕᴾ)(Bulk(2))),
        (ϕᶜ(Bulk(2)), bar(ϕᴾ)(Bulk(3))),
        (ϕᶜ(Bulk(3)), bar(ϕᴾ)(Bulk(1))),
    ]
    @test has_zero_loop(vs)

    vs = Contraction[
        (ϕᴾ(Bulk(1)), bar(ϕᶜ)(Bulk(2))),
        (ϕᴾ(Bulk(2)), bar(ϕᶜ)(Bulk(3))),
        (ϕᴾ(Bulk(3)), bar(ϕᶜ)(Bulk(1))),
    ]
    @test has_zero_loop(vs)

    vs = Contraction[
        (ϕᴾ(Bulk(1)), bar(ϕᶜ)(Bulk(2))),
        (ϕᶜ(Bulk(1)), bar(ϕᴾ)(Bulk(3))),
        (ϕᶜ(Bulk(3)), bar(ϕᴾ)(Bulk(2))),
    ]
    @test has_zero_loop(vs)

    vs = Contraction[
        (ϕᴾ(Bulk(1)), bar(ϕᶜ)(Bulk(2))),
        (ϕᶜ(Bulk(1)), bar(ϕᴾ)(Bulk(3))),
        (ϕᶜ(Bulk(2)), bar(ϕᴾ)(Bulk(3))),
    ]
    @test !has_zero_loop(vs)

    vs = Contraction[
        (ϕᶜ(Bulk(1)), bar(ϕᶜ)(Bulk(2))),
        (ϕᶜ(Bulk(2)), bar(ϕᶜ)(Bulk(3))),
        (ϕᶜ(Bulk(3)), bar(ϕᶜ)(Bulk(1))),
    ]
    @test !has_zero_loop(vs)

    retarded(i, j) = (ϕᶜ(Bulk(i)), bar(ϕᴾ)(Bulk(j)))
    advanced(i, j) = (ϕᴾ(Bulk(i)), bar(ϕᶜ)(Bulk(j)))
    keldysh(i, j) = (ϕᶜ(Bulk(i)), bar(ϕᶜ)(Bulk(j)))

    @testset "longer homogeneous cycles" begin
        @test has_zero_loop(
            Contraction[retarded(1, 2), retarded(2, 3), retarded(3, 4), retarded(4, 1)]
        )
        @test has_zero_loop(
            Contraction[
                advanced(1, 2),
                advanced(2, 3),
                advanced(3, 4),
                advanced(4, 5),
                advanced(5, 1),
            ],
        )
    end

    @testset "mixed cycles with extra edges" begin
        vs = Contraction[
            advanced(1, 2),
            retarded(3, 2),
            advanced(3, 4),
            retarded(5, 4),
            advanced(5, 1),
            advanced(1, 6),
        ]
        @test has_zero_loop(vs)

        vs = [vs[4], vs[2], vs[6], vs[1], vs[3], vs[5]]
        @test has_zero_loop(vs)
    end

    @testset "non-causal and acyclic graphs" begin
        @test !has_zero_loop(Contraction[advanced(1, 2), keldysh(2, 3), advanced(3, 1)])
        @test !has_zero_loop(
            Contraction[
                advanced(1, 2),
                retarded(3, 2),
                advanced(3, 4),
                retarded(5, 4),
                advanced(1, 5),
            ],
        )
    end
end
