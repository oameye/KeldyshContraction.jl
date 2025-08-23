using KeldyshContraction, Test
using KeldyshContraction: is_connected
using KeldyshContraction: is_irreducible, is_connected

@testset "Graph connectivity from multiplicity" begin
    @testset "Basic connectivity tests" begin
        # Single vertex
        @test is_connected(Int[], 1) == true

        # Two vertices, connected
        @test is_connected([1], 2) == true

        # Two vertices, disconnected
        @test is_connected([0], 2) == false

        # Triangle (complete graph K3)
        mult_triangle = [1, 1, 1]  # edges: (1,2), (1,3), (2,3)
        @test is_connected(mult_triangle, 3) == true

        # Path graph (1-2-3)
        mult_path = [1, 1, 0]  # edges: (1,2), (1,3) -> star with center 1
        @test is_connected(mult_path, 3) == true

        # Disconnected graph (two separate edges)
        mult_disconnected = [1, 0, 0, 0, 0, 1]  # 4 vertices: (1,2) and (3,4)
        @test is_connected(mult_disconnected, 4) == false
    end
end

@testset "Graph irreducibility" begin
    @testset "Basic irreducibility tests" begin
        # Single vertex - trivially irreducible
        @test is_irreducible(Int[], 1) == true

        # Two vertices - any single edge makes it reducible
        @test is_irreducible([1], 2) == false  # Bridge

        # Triangle (K3) - irreducible
        mult_triangle = [1, 1, 1]
        @test is_irreducible(mult_triangle, 3) == true

        # Path (1-2-3) - reducible (removing middle edge disconnects)
        mult_path = [1, 1, 0]
        @test is_irreducible(mult_path, 3) == false

        # Complete graph K4 - irreducible
        mult_k4 = [1, 1, 1, 1, 1, 1]
        @test is_irreducible(mult_k4, 4) == true

        # Square (cycle of 4) - irreducible
        mult_square = [1, 0, 1, 1, 0, 1]  # edges: (1,2), (1,4), (2,3), (3,4) -> cycle 1-2-3-4-1
        @test is_irreducible(mult_square, 4) == true
    end

    @testset "Edge cases" begin
        # Empty graph
        @test is_irreducible(Int[], 0) == true

        # Disconnected graph - should be irreducible (vacuously true)
        mult_disconnected = [1, 0, 0, 0, 0, 1]
        @test is_irreducible(mult_disconnected, 4) == false  # Not connected, so not irreducible

        # Multiple edges between same vertices
        mult_multi = [2, 1, 1]  # Double edge (1,2), others single
        @test is_irreducible(mult_multi, 3) == true  # Still irreducible
    end

    @testset "Bridge detection" begin
        # Star graph (one central vertex connected to others) - reducible
        mult_star = [1, 1, 1, 0, 0, 0]  # 4 vertices: 1 connected to 2,3,4
        @test is_irreducible(mult_star, 4) == false

        # Dumbbell graph (two triangles connected by bridge) - reducible
        # 5 vertices: triangle (1,2,3) + triangle (3,4,5)
        mult_dumbbell = [1, 1, 0, 0, 1, 1, 0, 0, 0, 1]  # Bridge at (3,4)
        @test is_irreducible(mult_dumbbell, 5) == false
    end
end

@testset "Performance and edge cases" begin
    @testset "Large graphs" begin
        # Complete graph K5
        mult_k5 = ones(Int, 10)  # 5 choose 2 = 10 edges
        @test is_irreducible(mult_k5, 5) == true

        # Path graph of length 5
        mult_path5 = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0]  # Only first 4 edges
        @test is_irreducible(mult_path5, 5) == false
    end

    @testset "Invalid inputs" begin
        # Empty multiplicity list for non-trivial graph should be handled gracefully
        @test is_connected(Int[], 2) == false
        @test is_irreducible(Int[], 2) == false
    end
end
