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
    # @test length(keys(topologies(GF3.keldysh))) == 10

    @test_broken GF = DressedPropagator(L_int, 4)
end

@testset "third order run's" begin
    @qfields c::Destroy(Classical) q::Destroy(Quantum)
    elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
    L_int = InteractionLagrangian(elasctic2boson)
    GF3 = DressedPropagator(L_int, 3)
    @test_broken SelfEnergy(GF3)
end

# using KeldyshContraction, Test

# using KeldyshContraction:
#     contractions,
#     is_irreducible,
#     position_category,
#     propagator_type,
#     is_retarded,
#     is_keldysh,
#     is_advanced,
#     SmallCollections,
#     PropagatorType,
#     is_bulk,
#     topology

# @qfields c::Destroy(Classical) q::Destroy(Quantum)
# elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
# L_int = InteractionLagrangian(elasctic2boson)
# GF3 = DressedPropagator(L_int, 3)

# CCcollections = Dict()
# E = typeof(GF3.keldysh).parameters[1]
# for (diagram, prefactor) in GF3.keldysh
#     _contractions = contractions(diagram)

#     if !is_irreducible(_contractions)
#         continue
#     end

#     positions = position_category.(_contractions)
#     types_p = propagator_type.(_contractions)
#     dict = SmallCollections.SmallDict{E,Symbol,PropagatorType.T}(
#         p => t for (p, t) in zip(positions, types_p)
#     )

#     # Find all bulk propagators (edges where both fields are bulk)
#     bulk_idxs = findall(is_bulk, _contractions)
#     bulk_propagators = _contractions[bulk_idxs]

#     if is_keldysh(dict[:out]) && is_advanced(dict[:in])
#         continue
#     elseif is_retarded(dict[:out]) && is_keldysh(dict[:in])
#         continue
#     elseif is_retarded(dict[:out]) && is_advanced(dict[:in])
#         continue
#     else
#         if haskey(CCcollections, diagram)
#             error("")
#         end
#         CCcollections[diagram] = prefactor
#     end
# end
# topology.(keys(CCcollections))

# topologies(GF3.keldysh)[[1, 1, 1]]
