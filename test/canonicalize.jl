using KeldyshContraction, Test
using KeldyshContraction: canonicalize, Bulk, In, Out, is_canonical

@qfields c::Destroy(Classical) q::Destroy(Quantum)

@testset "canonicalize" begin
    # Basic isomorphic graphs
    vs1 = [(c(Out()), q'(Bulk(1))), (c(Bulk(1)), q'(Bulk(2))), (c(Bulk(2)), q'(In()))]
    vs2 = [(c(Out()), q'(Bulk(2))), (c(Bulk(2)), q'(Bulk(1))), (c(Bulk(1)), q'(In()))]
    @test canonicalize(vs1) == canonicalize(vs2)

    # Original example - 3-node isomorphic graphs
    vs5 = [
        (c(Out()), q'(Bulk(1))),
        (c(Bulk(1)), q'(Bulk(2))),
        (c(Bulk(1)), q'(Bulk(3))),
        (c(Bulk(3)), q'(Bulk(1))),
        (c(Bulk(2)), q'(Bulk(3))),
        (c(Bulk(3)), q'(Bulk(2))),
        (c(Bulk(2)), q'(In())),
    ]
    vs6 = [
        (c(Out()), q'(Bulk(1))),
        (c(Bulk(1)), q'(Bulk(3))),
        (c(Bulk(1)), q'(Bulk(2))),
        (c(Bulk(2)), q'(Bulk(1))),
        (c(Bulk(3)), q'(Bulk(2))),
        (c(Bulk(2)), q'(Bulk(3))),
        (c(Bulk(3)), q'(In())),
    ]
    @test canonicalize(vs5) == canonicalize(vs6)

    # 4-node isomorphic graphs - ring topology with different node labeling
    vs_ring1 = [
        (c(Out()), q'(Bulk(1))),
        (c(Bulk(1)), q'(Bulk(2))),
        (c(Bulk(2)), q'(Bulk(3))),
        (c(Bulk(3)), q'(Bulk(4))),
        (c(Bulk(4)), q'(Bulk(1))), # ring closure
        (c(Bulk(4)), q'(In())),
    ]
    vs_ring2 = [
        (c(Out()), q'(Bulk(3))),
        (c(Bulk(3)), q'(Bulk(1))),
        (c(Bulk(1)), q'(Bulk(4))),
        (c(Bulk(4)), q'(Bulk(2))),
        (c(Bulk(2)), q'(Bulk(3))), # ring closure
        (c(Bulk(2)), q'(In())),
    ]
    @test canonicalize(vs_ring1) == canonicalize(vs_ring2)

    # Out() always connects to Bulk(1) in canonical form
    for vs in [vs1, vs2, vs5, vs6, vs_ring1, vs_ring2]
        canonical = canonicalize(vs)
        out_edge = findfirst(c -> Out() ∈ KeldyshContraction.position.(c), canonical)
        positions = KeldyshContraction.position.(canonical[out_edge])
        @test Bulk(1) ∈ positions
    end

    # Single bulk node (trivial case)
    vs_single = [(c(Out()), q'(Bulk(1))), (c(Bulk(1)), q'(In()))]
    @test canonicalize(vs_single) == vs_single

    # No bulk nodes (should return unchanged)
    vs_direct = [(c(Out()), q'(In()))]
    @test canonicalize(vs_direct) == vs_direct

    # Different graph topologies should have different canonical forms
    vs_linear = [(c(Out()), q'(Bulk(1))), (c(Bulk(1)), q'(Bulk(2))), (c(Bulk(2)), q'(In()))]
    vs_loop = [
        (c(Out()), q'(Bulk(1))),
        (c(Bulk(1)), q'(Bulk(2))),
        (c(Bulk(2)), q'(Bulk(1))), # loop between bulk nodes
        (c(Bulk(2)), q'(In())),
    ]
    @test canonicalize(vs_linear) != canonicalize(vs_loop)

    # Edge cases

    # Self-loops on bulk nodes
    vs_self_loop = [
        (c(Out()), q'(Bulk(1))),
        (c(Bulk(1)), q'(Bulk(1))), # self-loop
        (c(Bulk(1)), q'(In())),
    ]
    canonical_self = canonicalize(vs_self_loop)
    @test Bulk(1) ∈ KeldyshContraction.position.(canonical_self[1]) # Out() -> Bulk(1)

    # Multiple edges between same nodes
    vs_multi = [
        (c(Out()), q'(Bulk(1))),
        (c(Bulk(1)), q'(Bulk(2))),
        (c(Bulk(1)), q'(Bulk(2))), # duplicate edge
        (c(Bulk(2)), q'(In())),
    ]
    @test canonicalize(vs_multi) == canonicalize(vs_multi) # should be deterministic

    # Already canonical form (Out() connected to Bulk(1))
    vs_canonical = [
        (c(Out()), q'(Bulk(1))),
        (c(Bulk(1)), q'(Bulk(3))),
        (c(Bulk(3)), q'(Bulk(2))),
        (c(Bulk(2)), q'(In())),
    ]
    @test canonicalize(vs_canonical) == vs_canonical

    # Permutation invariance - same graph with different bulk numbering
    vs_perm1 = [
        (c(Out()), q'(Bulk(5))), (c(Bulk(5)), q'(Bulk(10))), (c(Bulk(10)), q'(In()))
    ]
    vs_perm2 = [
        (c(Out()), q'(Bulk(100))), (c(Bulk(100)), q'(Bulk(7))), (c(Bulk(7)), q'(In()))
    ]
    @test canonicalize(vs_perm1) == canonicalize(vs_perm2)
end
