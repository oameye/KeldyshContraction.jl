using KeldyshContraction, Test

@testset "number of topologies" begin
    @qfields c::Destroy(Classical) q::Destroy(Quantum)
    elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
    L_int = InteractionLagrangian(elasctic2boson)

    GF1 = DressedPropagator(L_int, 1)
    keys(topologies(GF1.keldysh))

    GF2 = DressedPropagator(L_int, 2)
    @test length(keys(topologies(GF2.keldysh))) == 3

    GF3 = DressedPropagator(L_int, 3)
    @test length(keys(topologies(GF3.keldysh))) == 10

    @test_broken GF = DressedPropagator(L_int, 4)
end
