using BenchmarkTools
using KeldyshContraction
using KeldyshContraction: canonicalize, Bulk, In, Out

function benchmark_canonicalize!(suite::BenchmarkGroup)
    @qfields c::Destroy(Classical) q::Destroy(Quantum)

    # 2-node linear graph
    vs_2node = [(c(Out()), q'(Bulk(1))), (c(Bulk(1)), q'(Bulk(2))), (c(Bulk(2)), q'(In()))]
    suite["microbenchmark"]["canonicalize"]["2-node"] = @benchmarkable canonicalize(
        $vs_2node
    )

    # 3-node fully connected graph
    vs_3node = [
        (c(Out()), q'(Bulk(2))),
        (c(Bulk(2)), q'(Bulk(1))),
        (c(Bulk(2)), q'(Bulk(3))),
        (c(Bulk(3)), q'(Bulk(2))),
        (c(Bulk(1)), q'(Bulk(3))),
        (c(Bulk(3)), q'(Bulk(1))),
        (c(Bulk(1)), q'(In())),
    ]
    suite["microbenchmark"]["canonicalize"]["3-node"] = @benchmarkable canonicalize(
        $vs_3node
    )

    # 4-node ring graph
    vs_4node_ring = [
        (c(Out()), q'(Bulk(3))),
        (c(Bulk(3)), q'(Bulk(1))),
        (c(Bulk(1)), q'(Bulk(4))),
        (c(Bulk(4)), q'(Bulk(2))),
        (c(Bulk(2)), q'(Bulk(3))),
        (c(Bulk(2)), q'(In())),
    ]
    suite["microbenchmark"]["canonicalize"]["4-node"] = @benchmarkable canonicalize(
        $vs_4node_ring
    )

    # 5-node complex graph (worst case - many permutations)
    vs_5node = [
        (c(Out()), q'(Bulk(5))),
        (c(Bulk(5)), q'(Bulk(1))),
        (c(Bulk(1)), q'(Bulk(2))),
        (c(Bulk(2)), q'(Bulk(3))),
        (c(Bulk(3)), q'(Bulk(4))),
        (c(Bulk(4)), q'(Bulk(5))),
        (c(Bulk(1)), q'(Bulk(3))),
        (c(Bulk(2)), q'(Bulk(4))),
        (c(Bulk(1)), q'(In())),
    ]
    suite["microbenchmark"]["canonicalize"]["5-node"] = @benchmarkable canonicalize(
        $vs_5node
    )

    return nothing
end
