using KeldyshContraction, Test

@testset "number of topologies" begin
    @qfields c::Destroy(Classical) q::Destroy(Quantum)
    elasctic2boson = -(1//2 * (c^2 + q^2) * c' * q' + 1//2 * c * q * ((c')^2 + (q')^2))
    L_int = InteractionLagrangian(elasctic2boson)

    GF1 = DressedPropagator(L_int, 1)
    @test_broken keys(topologies(GF1.keldysh))

    GF2 = DressedPropagator(L_int, 2)
    @test length(keys(topologies(GF2.keldysh))) == 3

    GF3 = DressedPropagator(L_int, 3)
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

    GF4 = DressedPropagator(L_int, 4);
    @test length(keys(topologies(GF4.keldysh))) == 59
    @test length(unique(sort.(keys(topologies(GF4.keldysh))))) == 17 # checked with mathematica (see test/All_graph_topologies.nb)

    irreduciable_topology = []
    for (key, value) in topologies(GF4.keldysh)
        cc = KeldyshContraction.contractions(first(value))
        if KeldyshContraction.is_irreducible(cc)
            push!(irreduciable_topology, key)
        end
    end
    length(unique(sort.(irreduciable_topology))) == 11 # checked with mathematica (see test/All_graph_topologies.nb)
end

@testset "third order run's" begin
    @qfields c::Destroy(Classical) q::Destroy(Quantum)
    elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
    L_int = InteractionLagrangian(elasctic2boson)
    GF3 = DressedPropagator(L_int, 3)
    @test_broken SelfEnergy(GF3)
end
