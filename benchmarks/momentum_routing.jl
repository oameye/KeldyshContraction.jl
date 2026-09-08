using KeldyshContraction
using KeldyshContraction: RoutingEdge, exact_momentum_routing, momentum_routing

function benchmark_momentum_routing!(SUITE)
    triangle = RoutingEdge[RoutingEdge(1, 2), RoutingEdge(2, 3), RoutingEdge(3, 1)]
    two_loop = RoutingEdge[
        RoutingEdge(1, 2),
        RoutingEdge(2, 3),
        RoutingEdge(3, 4),
        RoutingEdge(4, 1),
        RoutingEdge(1, 3),
    ]
    incidence = [
        -1 0 1 0 0
        1 -1 0 0 -1
        0 1 -1 -1 1
        0 0 0 1 0
    ]

    SUITE["Momentum routing"]["one loop graph"] = @benchmarkable momentum_routing($triangle) seconds =
        10
    SUITE["Momentum routing"]["two loop graph"] = @benchmarkable momentum_routing($two_loop) seconds =
        10
    SUITE["Momentum routing"]["exact matrix"] = @benchmarkable exact_momentum_routing(
        $incidence
    ) seconds = 10
    return nothing
end
