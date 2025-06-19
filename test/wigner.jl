using KeldyshContraction, Test
using KeldyshContraction:
    construct_linear_system, solve_linear_system, construct_momenta, momenta

@qfields c::Destroy(Classical) q::Destroy(Quantum)
elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
L_int = InteractionLagrangian(elasctic2boson)

@testset "Green's function" begin
    using KeldyshContraction: construct_momenta_from_gf
    @testset "Topology []" begin
        using KeldyshContraction: Momenta, FixedVector
        GF = DressedPropagator(L_int, 1)
        diagram = first(first(GF.keldysh))

        d′ = construct_momenta_from_gf(diagram)

        @test isequal(momenta(d′), FixedVector([Momenta(0), Momenta(1), Momenta(0)]))
    end

    GF = DressedPropagator(L_int, 2)

    topologies = KeldyshContraction.topologies(GF.keldysh)

    @testset "Topology [3]" begin
        diagram3 = first(topologies[[3]])
        A = construct_linear_system(diagram3.contractions)
        dep_idx, free_idx, P = solve_linear_system(A)
        @test dep_idx == [1, 4]
        @test free_idx == [2, 3, 5]
        @test P == [0.0 0.0 1.0; 1.0 1.0 -1.0]
        construct_momenta(dep_idx, free_idx, P)

        construct_momenta_from_gf(diagram3)
    end

    @testset "Topology [2]" begin
        diagram2 = first(topologies[[2]])
        A = construct_linear_system(diagram2.contractions)
        dep_idx, free_idx, P = solve_linear_system(A)
        @test dep_idx == [1, 3]
        @test free_idx == [2, 4, 5]
        @test P == [0.0 0.0 1.0; 1.0 0.0 0.0]
        construct_momenta(dep_idx, free_idx, P)

        construct_momenta_from_gf(diagram2)
    end

    @testset "Topology [1]" begin
        diagram1 = first(topologies[[1]])
        A = construct_linear_system(diagram1.contractions)
        dep_idx, free_idx, P = solve_linear_system(A)
        @test dep_idx == [1, 3]
        @test free_idx == [2, 4, 5]
        @test P == [0.0 0.0 1.0; 0.0 0.0 1.0]

        construct_momenta_from_gf(diagram1)
    end
end

@testset "Self-energy" begin
    using KeldyshContraction: construct_momenta_from_self_energy

    @testset "Topology []" begin
        using KeldyshContraction: Momenta, FixedVector
        GF = DressedPropagator(L_int, 1)
        SE = SelfEnergy(GF)

        diagram = first(first(SE.retarded))

        d′ = construct_momenta_from_self_energy(diagram)
        @test isequal(momenta(d′), FixedVector([Momenta(1)]))
    end

    GF = DressedPropagator(L_int, 2)
    SE = SelfEnergy(GF)

    topologies = KeldyshContraction.topologies(SE.retarded)

    @testset "Topology [3]" begin
        using KeldyshContraction: Momenta, FixedVector, Momentum
        diagram3 = first(topologies[[3]])
        d′ = construct_momenta_from_self_energy(diagram3)
        @test isequal(
            momenta(d′),
            FixedVector([
                Momenta(1),
                Momenta(2),
                Momenta([1, 1, -1], [Momentum(1), Momentum(2), Momentum(0)]),
            ]),
        )
    end

    @testset "Topology [2]" begin
        using KeldyshContraction: Momenta, FixedVector, Momentum
        diagram2 = first(topologies[[2]])
        d′ = construct_momenta_from_self_energy(diagram2)
        @test isequal(momenta(d′), FixedVector([
            Momenta(1),
            Momenta(1),
            Momenta(3), # rather have seen it to be Momenta(2)
        ]))
    end
end

# TODO topology -> positions
# function positions(topo::FixedVector{E,Int}) where {E}
#     return map(e -> e.position, topo)
# end
