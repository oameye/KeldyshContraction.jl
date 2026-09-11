using KeldyshContraction, Test
using KeldyshContraction:
    construct_linear_system, solve_linear_system, construct_momenta, momenta

@qfields ϕ::Boson
c, q = ϕ[Classical], ϕ[Quantum]
elasctic2boson = -(
    0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2)
)
L_int = InteractionLagrangian(elasctic2boson)

@testset "Green's function" begin
    using KeldyshContraction: construct_momenta_from_gf
    @testset "Topology []" begin
        using KeldyshContraction: Momenta, FixedVector
        GF = DressedPropagator(L_int, Val(1), Val(3))
        diagram = first(first(GF.keldysh))

        d′ = @inferred construct_momenta_from_gf(diagram)

        @test isequal(momenta(d′), FixedVector([Momenta(0), Momenta(1), Momenta(0)]))
    end

    GF = DressedPropagator(L_int, Val(2), Val(5))

    topologies = KeldyshContraction.topologies(GF.keldysh)

    @testset "Topology [3]" begin
        diagram3 = first(topologies[[3]])
        A = @inferred construct_linear_system(diagram3.contractions)
        dep_idx, free_idx, P = @inferred solve_linear_system(A)
        @test dep_idx == [1, 4]
        @test free_idx == [2, 3, 5]
        @test P == [0.0 0.0 1.0; 1.0 1.0 -1.0]
        @inferred construct_momenta(dep_idx, free_idx, P)

        @inferred construct_momenta_from_gf(diagram3)
    end

    @testset "Topology [2]" begin
        diagram2 = first(topologies[[2]])
        A = @inferred construct_linear_system(diagram2.contractions)
        dep_idx, free_idx, P = @inferred solve_linear_system(A)
        @test dep_idx == [1, 3]
        @test free_idx == [2, 4, 5]
        @test P == [0.0 0.0 1.0; 1.0 0.0 0.0]
        @inferred construct_momenta(dep_idx, free_idx, P)

        @inferred construct_momenta_from_gf(diagram2)
    end

    @testset "Topology [1]" begin
        diagram1 = first(topologies[[1]])
        A = @inferred construct_linear_system(diagram1.contractions)
        dep_idx, free_idx, P = @inferred solve_linear_system(A)
        @test dep_idx == [1, 3]
        @test free_idx == [2, 4, 5]
        @test P == [0.0 0.0 1.0; 0.0 0.0 1.0]

        @inferred construct_momenta_from_gf(diagram1)
    end
end

@testset "Self-energy" begin
    using KeldyshContraction: construct_momenta_from_self_energy

    @testset "Topology []" begin
        using KeldyshContraction: Momenta, FixedVector
        GF = DressedPropagator(L_int, Val(1), Val(3))
        SE = @inferred SelfEnergy(GF)
        @test SE isa SelfEnergy{ComplexF64,Boson,1,1,0}

        diagram = first(first(SE.retarded))

        d′ = @inferred construct_momenta_from_self_energy(diagram)
        @test isequal(momenta(d′), FixedVector([Momenta(1)]))
    end

    GF = DressedPropagator(L_int, Val(2), Val(5))
    SE = SelfEnergy(GF)

    topologies = KeldyshContraction.topologies(SE.retarded)

    @testset "Topology [3]" begin
        using KeldyshContraction: Momenta, FixedVector, Momentum
        diagram3 = first(topologies[[3]])
        d′ = @inferred construct_momenta_from_self_energy(diagram3)
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
        d′ = @inferred construct_momenta_from_self_energy(diagram2)
        @test isequal(momenta(d′), FixedVector([Momenta(1), Momenta(1), Momenta(3)]))
    end
end

@testset "Wigner transform" begin
    @testset "first order" begin
        GF = DressedPropagator(L_int, Val(1), Val(3))
        SE = SelfEnergy(GF)

        GFk = @inferred wigner_transform(GF)
        SEk = @inferred wigner_transform(SE)
        @test typeof(GFk) === typeof(GF)
        @test typeof(SEk) === typeof(SE)
    end
    @testset "second order" begin
        GF = DressedPropagator(L_int, Val(2), Val(5))
        SE = SelfEnergy(GF)

        GFk = @inferred wigner_transform(GF)
        SEk = @inferred wigner_transform(SE)
        @test typeof(GFk) === typeof(GF)
        @test typeof(SEk) === typeof(SE)
    end
end
