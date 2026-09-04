using KeldyshContraction, Test

function has_no_zero_loops(component)
    return all(component) do pair
        diagram = first(pair)
        vs = KeldyshContraction.Contraction[
            (edge.out, edge.in) for edge in KeldyshContraction.contractions(diagram)
        ]
        return !KeldyshContraction.has_zero_loop(vs)
    end
end

function has_valid_self_energy_diagrams(component)
    return all(component) do pair
        edges = KeldyshContraction.contractions(first(pair))
        return all(KeldyshContraction.is_bulk, edges) &&
               KeldyshContraction.is_irreducible(edges)
    end
end

@testset "number of topologies" begin
    @qfields c::Boson(Classical) q::Boson(Quantum)
    elasctic2boson = -(
        1//2 * (c^2 + q^2) * bar(c) * bar(q) +
        1//2 * c * q * (bar(c)^2 + bar(q)^2)
    )
    L_int = InteractionLagrangian(elasctic2boson)

    GF1 = DressedPropagator(L_int, Val(1), Val(3))
    @test length(keys(topologies(GF1.keldysh))) == 1

    GF2 = DressedPropagator(L_int, Val(2), Val(5))
    @test length(keys(topologies(GF2.keldysh))) == 3

    GF3 = DressedPropagator(L_int, Val(3), Val(7))
    @test length(keys(topologies(GF3.keldysh))) == 11
    @test length(unique(sort.(keys(topologies(GF3.keldysh))))) == 8

    irreduciable_topology = []
    for (key, value) in topologies(GF3.keldysh)
        cc = KeldyshContraction.contractions(first(value))
        if KeldyshContraction.is_irreducible(cc)
            push!(irreduciable_topology, key)
        end
    end
    @test length(unique(sort.(irreduciable_topology))) == 5

    @testset "zero-loop filtering" begin
        @test all(has_no_zero_loops, (GF3.keldysh, GF3.retarded, GF3.advanced))
    end

    GF4 = DressedPropagator(L_int, Val(4), Val(9))
    @test length(keys(topologies(GF4.keldysh))) == 59
    @test length(unique(sort.(keys(topologies(GF4.keldysh))))) == 17

    irreduciable_topology = []
    for (key, value) in topologies(GF4.keldysh)
        cc = KeldyshContraction.contractions(first(value))
        if KeldyshContraction.is_irreducible(cc)
            push!(irreduciable_topology, key)
        end
    end
    @test length(unique(sort.(irreduciable_topology))) == 11

    @testset "zero-loop filtering" begin
        @test all(has_no_zero_loops, (GF4.keldysh, GF4.retarded, GF4.advanced))
    end
end

@testset "third order run's" begin
    @qfields c::Boson(Classical) q::Boson(Quantum)
    elasctic2boson = -(
        0.5 * (c^2 + q^2) * bar(c) * bar(q) +
        0.5 * c * q * (bar(c)^2 + bar(q)^2)
    )
    L_int = InteractionLagrangian(elasctic2boson)
    GF3 = DressedPropagator(L_int, Val(3), Val(7))
    Σ = SelfEnergy(GF3, Val(3))
    @test Σ isa SelfEnergy
    @test all(component -> !isempty(component), (Σ.keldysh, Σ.retarded, Σ.advanced))
    @test all(has_valid_self_energy_diagrams, (Σ.keldysh, Σ.retarded, Σ.advanced))
end
