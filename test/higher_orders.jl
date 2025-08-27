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

    GF4 = DressedPropagator(L_int, 4)
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

using KeldyshContraction, Test

using KeldyshContraction:
    contractions,
    is_irreducible,
    position_category,
    propagator_type,
    is_retarded,
    is_keldysh,
    is_advanced,
    SmallCollections,
    PropagatorType,
    is_bulk,
    topology

@qfields c::Destroy(Classical) q::Destroy(Quantum)
elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
L_int = InteractionLagrangian(elasctic2boson)
GF3 = DressedPropagator(L_int, 3)

CCcollections = Dict()
E = typeof(GF3.keldysh).parameters[1]
for (diagram, prefactor) in GF3.keldysh
    _contractions = contractions(diagram)

    if !is_irreducible(_contractions)
        continue
    end

    positions = position_category.(_contractions)
    types_p = propagator_type.(_contractions)
    dict = SmallCollections.SmallDict{E,Symbol,PropagatorType.T}(
        p => t for (p, t) in zip(positions, types_p)
    )

    # Find all bulk propagators (edges where both fields are bulk)
    bulk_idxs = findall(is_bulk, _contractions)
    bulk_propagators = _contractions[bulk_idxs]

    if is_keldysh(dict[:out]) && is_advanced(dict[:in])
        continue
    elseif is_retarded(dict[:out]) && is_keldysh(dict[:in])
        continue
    elseif is_retarded(dict[:out]) && is_advanced(dict[:in])
        continue
    else
        if haskey(CCcollections, diagram)
            error("")
        end
        CCcollections[diagram] = prefactor
    end
end
CCcollections
unique(topology.(keys(CCcollections)))

idxs = findall(pair -> topology(pair[1]) == [2, 1, 1], collect(CCcollections))
collect(CCcollections)[idxs] # are zero due loop with only G^A or G^R


idxs = findall(pair -> topology(pair[1]) == [1, 3, 1], collect(CCcollections))
collect(CCcollections)[idxs] # are zero due loop with only G^A or G^R


test = (x -> (x.out, x.in)).(contractions(collect(CCcollections)[idxs][2][1]))
test_v = [x for x in test]

_has_zero_loop(test_v)

g, max_label, has_in = make_graph(Graphs.SimpleDiGraph, test_v)
simplecycles_iter(g)
