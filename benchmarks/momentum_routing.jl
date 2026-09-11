using KeldyshContraction
using KeldyshContraction:
    Bulk,
    Contraction,
    Diagram,
    FourierDiagram,
    In,
    Out,
    RoutingEdge,
    exact_affine_momentum_routing,
    exact_momentum_routing,
    lower_fourier_derivatives,
    momentum_routing

@qfields benchmark_fourier_ϕ::Boson
const benchmark_fourier_c = benchmark_fourier_ϕ[Classical]
const benchmark_fourier_q = benchmark_fourier_ϕ[Quantum]

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
    affine_incidence = [-1 0; 1 -1; 0 1]
    affine_source = reshape([-1, 0, 1], 3, 1)

    contractions = Contraction{Boson}[
        Contraction(benchmark_fourier_c(Out()), bar(benchmark_fourier_q)(Bulk(1))),
        Contraction(benchmark_fourier_c(Bulk(2)), bar(benchmark_fourier_q)(Bulk(1))),
        Contraction(benchmark_fourier_c(Bulk(2)), bar(benchmark_fourier_q)(Bulk(1))),
        Contraction(benchmark_fourier_c(Bulk(1)), bar(benchmark_fourier_q)(Bulk(2))),
        Contraction(benchmark_fourier_c(Bulk(2)), bar(benchmark_fourier_q)(In())),
    ]
    diagram = Diagram(contractions, Val(5), Val(1))

    derivative_contractions = Contraction{Boson}[
        Contraction(
            partial(benchmark_fourier_c, :x)(Out()), bar(benchmark_fourier_q)(Bulk(1))
        ),
        Contraction(benchmark_fourier_c(Bulk(1)), bar(benchmark_fourier_q)(Bulk(1))),
        Contraction(benchmark_fourier_c(Bulk(1)), bar(benchmark_fourier_q)(In())),
    ]
    derivative_diagram = Diagram(derivative_contractions, Val(3), Val(0))
    routed_derivative = FourierDiagram(derivative_diagram)

    SUITE["Momentum routing"]["one loop graph"] = @benchmarkable momentum_routing($triangle) seconds =
        10
    SUITE["Momentum routing"]["two loop graph"] = @benchmarkable momentum_routing($two_loop) seconds =
        10
    SUITE["Momentum routing"]["exact matrix"] = @benchmarkable exact_momentum_routing(
        $incidence
    ) seconds = 10
    SUITE["Momentum routing"]["exact affine matrix"] = @benchmarkable exact_affine_momentum_routing(
        $affine_incidence, $affine_source
    ) seconds = 10
    SUITE["Momentum routing"]["Fourier diagram"] = @benchmarkable FourierDiagram($diagram) seconds =
        10
    SUITE["Momentum routing"]["derivative lowering"] = @benchmarkable lower_fourier_derivatives(
        $routed_derivative
    ) seconds = 10
    return nothing
end
