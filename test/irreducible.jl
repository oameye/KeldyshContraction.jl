using KeldyshContraction, Test
using KeldyshContraction: is_connected, is_irreducible

@testset "Irreducibility for Vector{Contraction}" begin
    using KeldyshContraction: Contraction, In, Out, Bulk

    @qfields c::Boson(Classical) q::Boson(Quantum)

    @testset "Basic contraction irreducibility" begin
        @test is_irreducible(Contraction[]) == true

        vs_single = Contraction[(c, bar(q))]
        @test is_irreducible(vs_single) == true

        vs_connected = Contraction[(c(Out()), bar(q)), (c, bar(q)), (c, bar(q)(In()))]
        @test is_connected(vs_connected)
        @test is_irreducible(vs_connected) == true
    end

    @testset "Bridge detection with multiple vertices" begin
        vs_dumbbell = Contraction[
            (c(Bulk(1)), bar(q)(Bulk(2))),
            (c(Bulk(2)), bar(q)(Bulk(3))),
            (c(Bulk(3)), bar(q)(Bulk(1))),
            (c(Bulk(3)), bar(q)(Bulk(4))),
            (c(Bulk(4)), bar(q)(Bulk(5))),
            (c(Bulk(5)), bar(q)(Bulk(6))),
            (c(Bulk(6)), bar(q)(Bulk(4))),
        ]

        @test is_connected(vs_dumbbell)
        @test !is_irreducible(vs_dumbbell)
    end

    @testset "Triangle graph - irreducible" begin
        vs_triangle = Contraction[
            (c(Bulk(1)), bar(q)(Bulk(2))),
            (c(Bulk(2)), bar(q)(Bulk(3))),
            (c(Bulk(3)), bar(q)(Bulk(1))),
        ]

        @test is_connected(vs_triangle)
        @test is_irreducible(vs_triangle)
    end

    @testset "Square with diagonal - irreducible" begin
        vs_square_diag = Contraction[
            (c(Bulk(1)), bar(q)(Bulk(2))),
            (c(Bulk(2)), bar(q)(Bulk(3))),
            (c(Bulk(3)), bar(q)(Bulk(4))),
            (c(Bulk(4)), bar(q)(Bulk(1))),
            (c(Bulk(1)), bar(q)(Bulk(3))),
        ]

        @test is_connected(vs_square_diag)
        @test is_irreducible(vs_square_diag)
    end

    @testset "Complete graph K4 - irreducible" begin
        vs_k4 = Contraction[
            (c(Bulk(1)), bar(q)(Bulk(2))),
            (c(Bulk(1)), bar(q)(Bulk(3))),
            (c(Bulk(1)), bar(q)(Bulk(4))),
            (c(Bulk(2)), bar(q)(Bulk(3))),
            (c(Bulk(2)), bar(q)(Bulk(4))),
            (c(Bulk(3)), bar(q)(Bulk(4))),
        ]

        @test is_connected(vs_k4)
        @test is_irreducible(vs_k4)
    end

    @testset "Edge cases" begin
        vs_self = Contraction[(c(Out()), bar(c)(In()))]
        @test is_irreducible(vs_self) == true
        @test is_irreducible(Contraction[]) == true
    end
end
