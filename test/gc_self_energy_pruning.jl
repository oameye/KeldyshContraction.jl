using KeldyshContraction: _has_unhealable_bulk_bridge
using Test

@testset "Continuation-safe 1PI bridge pruning" begin
    # A realized chain has an unavoidable bridge when no residual edge can heal either cut.
    @test _has_unhealable_bulk_bridge([(1, 2), (2, 3)], Tuple{Int,Int}[])

    # Current reducibility is not monotone: the residual 1--3 edge heals both realized bridges.
    @test !_has_unhealable_bulk_bridge([(1, 2), (2, 3)], [(1, 3)])

    # Already bridgeless partial support must never be rejected.
    @test !_has_unhealable_bulk_bridge([(1, 2), (2, 3), (3, 1)], Tuple{Int,Int}[])

    # A parallel realized propagator heals the other copy exactly as in the final multigraph.
    @test !_has_unhealable_bulk_bridge([(1, 2), (1, 2)], Tuple{Int,Int}[])

    # An unrelated optional edge cannot heal a realized bridge.
    @test _has_unhealable_bulk_bridge([(1, 2), (2, 3)], [(3, 4)])

    # Optional dangling support must not make an otherwise bridgeless partial graph rejectable.
    @test !_has_unhealable_bulk_bridge([(1, 2), (2, 3), (3, 1)], [(3, 4)])

    # Self-loops are never bridges under the existing KC irreducibility convention.
    @test !_has_unhealable_bulk_bridge([(1, 1)], Tuple{Int,Int}[])
end
