using KeldyshContraction, Test
using KeldyshContraction: canonicalize, Bulk, In, Out, sort_by_position_and_type, positions

@qfields ϕ::Boson
c, q = ϕ[Classical], ϕ[Quantum]

@testset "canonicalize" begin
    # Basic isomorphic graphs
    vs1 = [
        (c(Out()), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(Bulk(2))),
        (c(Bulk(2)), bar(q)(In())),
    ]
    vs2 = [
        (c(Out()), bar(q)(Bulk(2))),
        (c(Bulk(2)), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(In())),
    ]
    @test canonicalize(vs1) == canonicalize(vs2)

    # Original example - 3-node isomorphic graphs
    vs5 = [
        (c(Out()), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(Bulk(2))),
        (c(Bulk(1)), bar(q)(Bulk(3))),
        (c(Bulk(3)), bar(q)(Bulk(1))),
        (c(Bulk(2)), bar(q)(Bulk(3))),
        (c(Bulk(3)), bar(q)(Bulk(2))),
        (c(Bulk(2)), bar(q)(In())),
    ]
    vs6 = [
        (c(Out()), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(Bulk(3))),
        (c(Bulk(1)), bar(q)(Bulk(2))),
        (c(Bulk(2)), bar(q)(Bulk(1))),
        (c(Bulk(3)), bar(q)(Bulk(2))),
        (c(Bulk(2)), bar(q)(Bulk(3))),
        (c(Bulk(3)), bar(q)(In())),
    ]
    @test canonicalize(vs5) == canonicalize(vs6)

    # 4-node isomorphic graphs - ring topology with different node labeling
    vs_ring1 = [
        (c(Out()), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(Bulk(2))),
        (c(Bulk(2)), bar(q)(Bulk(3))),
        (c(Bulk(3)), bar(q)(Bulk(4))),
        (c(Bulk(4)), bar(q)(Bulk(1))),
        (c(Bulk(4)), bar(q)(In())),
    ]
    vs_ring2 = [
        (c(Out()), bar(q)(Bulk(3))),
        (c(Bulk(3)), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(Bulk(4))),
        (c(Bulk(4)), bar(q)(Bulk(2))),
        (c(Bulk(2)), bar(q)(Bulk(3))),
        (c(Bulk(2)), bar(q)(In())),
    ]
    @test canonicalize(vs_ring1) == canonicalize(vs_ring2)

    # Out() always connects to Bulk(1) in canonical form
    for vs in [vs1, vs2, vs5, vs6, vs_ring1, vs_ring2]
        canonical = canonicalize(vs)
        out_edge = findfirst(cn -> Out() ∈ KeldyshContraction.position.(cn), canonical)
        ps = KeldyshContraction.position.(canonical[out_edge])
        @test Bulk(1) ∈ ps
    end

    vs_single = [(c(Out()), bar(q)(Bulk(1))), (c(Bulk(1)), bar(q)(In()))]
    @test canonicalize(vs_single) == vs_single

    vs_direct = [(c(Out()), bar(q)(In()))]
    @test canonicalize(vs_direct) == vs_direct

    vs_linear = [
        (c(Out()), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(Bulk(2))),
        (c(Bulk(2)), bar(q)(In())),
    ]
    vs_loop = [
        (c(Out()), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(Bulk(2))),
        (c(Bulk(2)), bar(q)(Bulk(1))),
        (c(Bulk(2)), bar(q)(In())),
    ]
    @test canonicalize(vs_linear) != canonicalize(vs_loop)

    vs_self_loop = [
        (c(Out()), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(In())),
    ]
    canonical_self = canonicalize(vs_self_loop)
    @test Bulk(1) ∈ KeldyshContraction.position.(canonical_self[1])

    vs_multi = [
        (c(Out()), bar(q)(Bulk(1))),
        (c(Bulk(1)), bar(q)(Bulk(2))),
        (c(Bulk(1)), bar(q)(Bulk(2))),
        (c(Bulk(2)), bar(q)(In())),
    ]
    @test canonicalize(vs_multi) == canonicalize(vs_multi)

    @testset "third order two body scattering" begin
        vs1 = [
            (c(Out()), bar(q)(Bulk(1))),
            (c(Bulk(1)), bar(q)(Bulk(2))),
            (c(Bulk(2)), bar(q)(Bulk(2))),
            (c(Bulk(3)), bar(q)(Bulk(1))),
            (c(Bulk(2)), bar(q)(Bulk(3))),
            (c(Bulk(3)), bar(q)(Bulk(3))),
            (c(Bulk(1)), bar(q)(In())),
        ]
        vs2 = [
            (c(Out()), bar(q)(Bulk(2))),
            (c(Bulk(1)), bar(q)(Bulk(1))),
            (c(Bulk(1)), bar(q)(Bulk(2))),
            (c(Bulk(3)), bar(q)(Bulk(1))),
            (c(Bulk(2)), bar(q)(Bulk(3))),
            (c(Bulk(3)), bar(q)(Bulk(3))),
            (c(Bulk(2)), bar(q)(In())),
        ]
        vscanonical1 = canonicalize(vs1)
        vscanonical2 = canonicalize(vs2)
        sort!(vscanonical1; by=sort_by_position_and_type)
        sort!(vscanonical2; by=sort_by_position_and_type)
        @test vscanonical1 == vscanonical2
    end

    @testset "type_stability" begin
        using KeldyshContraction: make_NautyDiGraph

        vs1 = [
            (c(Out()), bar(q)(Bulk(1))),
            (c(Bulk(1)), bar(q)(Bulk(2))),
            (c(Bulk(2)), bar(q)(Bulk(2))),
            (c(Bulk(3)), bar(q)(Bulk(1))),
            (c(Bulk(2)), bar(q)(Bulk(3))),
            (c(Bulk(3)), bar(q)(Bulk(3))),
            (c(Bulk(1)), bar(q)(In())),
        ]
        @inferred make_NautyDiGraph(vs1)
    end
end
