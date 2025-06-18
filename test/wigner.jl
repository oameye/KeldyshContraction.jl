using KeldyshContraction, Test
using KeldyshContraction: construct_linear_system, solve_linear_system

@qfields c::Destroy(Classical) q::Destroy(Quantum)
elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
L_int = InteractionLagrangian(elasctic2boson)

GF = DressedPropagator(L_int, 2)

topologies = KeldyshContraction.topologies(GF.keldysh)

@testset "Topology [3]" begin
    diagram3 = first(topologies[[3]])
    A = construct_linear_system(diagram3.contractions)
    dep_idx, free_idx, P = solve_linear_system(A)
    @test dep_idx == [1, 4]
    @test free_idx == [2, 3, 5]
    @test P == [0.0 0.0 1.0; 1.0 1.0 -1.0]
end

@testset "Topology [2]" begin
    diagram2 = first(topologies[[2]])
    A = construct_linear_system(diagram2.contractions)
    dep_idx, free_idx, P = solve_linear_system(A)
    @test dep_idx == [1, 3]
    @test free_idx == [2, 4, 5]
    @test P == [0.0 0.0 1.0; 1.0 0.0 0.0]
end

@testset "Topology [1]" begin
    diagram1 = first(topologies[[1]])
    A = construct_linear_system(diagram1.contractions)
    dep_idx, free_idx, P = solve_linear_system(A)
    @test dep_idx == [1, 3]
    @test free_idx == [2, 4, 5]
    @test P == [0.0 0.0 1.0; 0.0 0.0 1.0]
end
