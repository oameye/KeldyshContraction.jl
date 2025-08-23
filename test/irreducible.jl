using KeldyshContraction, Test
using KeldyshContraction: is_connected, is_irreducible

@testset "Irreducibility for Vector{Contraction}" begin
    using KeldyshContraction: Contraction, In, Out, Bulk

    @qfields c::Destroy(Classical) q::Destroy(Quantum)

    @testset "Basic contraction irreducibility" begin
        # Empty vector - trivially irreducible
        @test is_irreducible(Contraction[]) == true

        # Single contraction - trivially irreducible
        vs_single = Contraction[(c, q')]
        @test is_irreducible(vs_single) == true

        # Known connected example from diagram.jl tests
        vs_connected = Contraction[(c(Out()), q'), (c, q'), (c, q(In())')]
        @test is_connected(vs_connected)  # Verify this is connected
        @test is_irreducible(vs_connected) == true  # Should be irreducible (triangle-like)
    end

    @testset "Bridge detection with multiple vertices" begin
        # Create a "dumbbell" graph: triangle-bridge-triangle
        # Vertices: 1-2-3 form triangle, 3-4 is bridge, 4-5-6 form triangle
        # This should be reducible (removing the bridge 3-4 disconnects the triangles)
        vs_dumbbell = Contraction[
            (c(Bulk(1)), q'(Bulk(2))),  # Triangle 1: edge 1-2
            (c(Bulk(2)), q'(Bulk(3))),  # Triangle 1: edge 2-3
            (c(Bulk(3)), q'(Bulk(1))),  # Triangle 1: edge 3-1
            (c(Bulk(3)), q'(Bulk(4))),  # Bridge: edge 3-4
            (c(Bulk(4)), q'(Bulk(5))),  # Triangle 2: edge 4-5
            (c(Bulk(5)), q'(Bulk(6))),  # Triangle 2: edge 5-6
            (c(Bulk(6)), q'(Bulk(4))),   # Triangle 2: edge 6-4
        ]

        @test is_connected(vs_dumbbell)  # Should be connected
        @test !is_irreducible(vs_dumbbell)  # Should be reducible (removing bridge disconnects)
    end

    @testset "Triangle graph - irreducible" begin
        # Create a triangle: vertex 1 - vertex 2 - vertex 3 - vertex 1
        # This should be irreducible (no single edge removal disconnects it)
        vs_triangle = Contraction[
            (c(Bulk(1)), q'(Bulk(2))),  # Edge 1-2
            (c(Bulk(2)), q'(Bulk(3))),  # Edge 2-3
            (c(Bulk(3)), q'(Bulk(1))),   # Edge 3-1
        ]

        @test is_connected(vs_triangle)  # Should be connected
        @test is_irreducible(vs_triangle)  # Should be irreducible (triangle)
    end

    @testset "Square with diagonal - irreducible" begin
        # Create a square with a diagonal: 1-2-3-4-1 plus diagonal 1-3
        # This should be irreducible (removing any single edge still leaves connectivity)
        vs_square_diag = Contraction[
            (c(Bulk(1)), q'(Bulk(2))),  # Square: edge 1-2
            (c(Bulk(2)), q'(Bulk(3))),  # Square: edge 2-3
            (c(Bulk(3)), q'(Bulk(4))),  # Square: edge 3-4
            (c(Bulk(4)), q'(Bulk(1))),  # Square: edge 4-1
            (c(Bulk(1)), q'(Bulk(3))),   # Diagonal: edge 1-3
        ]

        @test is_connected(vs_square_diag)  # Should be connected
        @test is_irreducible(vs_square_diag)  # Should be irreducible (multiple paths between vertices)
    end

    @testset "Complete graph K4 - irreducible" begin
        # Create a complete graph on 4 vertices
        # This should be irreducible
        vs_k4 = Contraction[
            (c(Bulk(1)), q'(Bulk(2))),  # Edge 1-2
            (c(Bulk(1)), q'(Bulk(3))),  # Edge 1-3
            (c(Bulk(1)), q'(Bulk(4))),  # Edge 1-4
            (c(Bulk(2)), q'(Bulk(3))),  # Edge 2-3
            (c(Bulk(2)), q'(Bulk(4))),  # Edge 2-4
            (c(Bulk(3)), q'(Bulk(4))),   # Edge 3-4
        ]

        @test is_connected(vs_k4)  # Should be connected
        @test is_irreducible(vs_k4)  # Should be irreducible (complete graphs are irreducible)
    end

    @testset "Edge cases" begin
        # Single self-loop type contraction
        vs_self = Contraction[(c(Out()), c'(In()))]
        @test is_irreducible(vs_self) == true  # Single contraction is trivially irreducible

        # Test the empty case explicitly
        @test is_irreducible(Contraction[]) == true
    end
end
