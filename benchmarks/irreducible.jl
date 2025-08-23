using BenchmarkTools
using KeldyshContraction
using KeldyshContraction: is_irreducible, Bulk, In, Out

function benchmark_irreducible!(suite::BenchmarkGroup)
    @qfields c::Destroy(Classical) q::Destroy(Quantum)

    # Dumbbell graph (6 vertices) - reducible with bridge
    vs_dumbbell = [
        (c(Bulk(1)), q'(Bulk(2))),  # Triangle 1: edge 1-2
        (c(Bulk(2)), q'(Bulk(3))),  # Triangle 1: edge 2-3
        (c(Bulk(3)), q'(Bulk(1))),  # Triangle 1: edge 3-1
        (c(Bulk(3)), q'(Bulk(4))),  # Bridge: edge 3-4
        (c(Bulk(4)), q'(Bulk(5))),  # Triangle 2: edge 4-5
        (c(Bulk(5)), q'(Bulk(6))),  # Triangle 2: edge 5-6
        (c(Bulk(6)), q'(Bulk(4))),   # Triangle 2: edge 6-4
    ]
    suite["microbenchmark"]["irreducible"]["dumbbell"] = @benchmarkable is_irreducible(
        $vs_dumbbell
    )

    # Complete graph K4 (4 vertices) - irreducible
    vs_k4 = [
        (c(Bulk(1)), q'(Bulk(2))),  # Edge 1-2
        (c(Bulk(1)), q'(Bulk(3))),  # Edge 1-3
        (c(Bulk(1)), q'(Bulk(4))),  # Edge 1-4
        (c(Bulk(2)), q'(Bulk(3))),  # Edge 2-3
        (c(Bulk(2)), q'(Bulk(4))),  # Edge 2-4
        (c(Bulk(3)), q'(Bulk(4))),   # Edge 3-4
    ]
    suite["microbenchmark"]["irreducible"]["complete-k4"] = @benchmarkable is_irreducible(
        $vs_k4
    )

    # Square with diagonal (4 vertices) - irreducible
    vs_square_diag = [
        (c(Bulk(1)), q'(Bulk(2))),  # Square: edge 1-2
        (c(Bulk(2)), q'(Bulk(3))),  # Square: edge 2-3
        (c(Bulk(3)), q'(Bulk(4))),  # Square: edge 3-4
        (c(Bulk(4)), q'(Bulk(1))),  # Square: edge 4-1
        (c(Bulk(1)), q'(Bulk(3))),   # Diagonal: edge 1-3
    ]
    suite["microbenchmark"]["irreducible"]["square-diagonal"] = @benchmarkable is_irreducible(
        $vs_square_diag
    )

    # Complete graph K5 (5 vertices) - irreducible, more expensive
    vs_k5 = [
        (c(Bulk(1)), q'(Bulk(2))),  # Edge 1-2
        (c(Bulk(1)), q'(Bulk(3))),  # Edge 1-3
        (c(Bulk(1)), q'(Bulk(4))),  # Edge 1-4
        (c(Bulk(1)), q'(Bulk(5))),  # Edge 1-5
        (c(Bulk(2)), q'(Bulk(3))),  # Edge 2-3
        (c(Bulk(2)), q'(Bulk(4))),  # Edge 2-4
        (c(Bulk(2)), q'(Bulk(5))),  # Edge 2-5
        (c(Bulk(3)), q'(Bulk(4))),  # Edge 3-4
        (c(Bulk(3)), q'(Bulk(5))),  # Edge 3-5
        (c(Bulk(4)), q'(Bulk(5))),   # Edge 4-5
    ]
    suite["microbenchmark"]["irreducible"]["complete-k5"] = @benchmarkable is_irreducible(
        $vs_k5
    )

    # Large sparse reducible graph (7 vertices) - tree structure
    vs_tree = [
        (c(Bulk(1)), q'(Bulk(2))),  # Root to left subtree
        (c(Bulk(1)), q'(Bulk(3))),  # Root to right subtree
        (c(Bulk(2)), q'(Bulk(4))),  # Left subtree
        (c(Bulk(2)), q'(Bulk(5))),  # Left subtree
        (c(Bulk(3)), q'(Bulk(6))),  # Right subtree
        (c(Bulk(3)), q'(Bulk(7))),   # Right subtree
    ]
    suite["microbenchmark"]["irreducible"]["tree"] = @benchmarkable is_irreducible($vs_tree)

    return nothing
end
