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

    # 4-node ring regression for #184: bulk labels must canonicalize completely.
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

    # External labels are fixed: Out() remains attached to canonical Bulk(1).
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

    @testset "color-aware bulk relabeling" begin
        using KeldyshContraction: FieldIndex, FieldIndices
        @qfields χ::Boson
        χc, χq = χ[Classical], χ[Quantum]

        colored1 = [
            (c(Out()), bar(q)(Bulk(1))),
            (χc(Bulk(1)), bar(χq)(Bulk(2))),
            (c(Bulk(2)), bar(q)(In())),
        ]
        colored2 = [
            (c(Out()), bar(q)(Bulk(2))),
            (χc(Bulk(2)), bar(χq)(Bulk(1))),
            (c(Bulk(1)), bar(q)(In())),
        ]
        @test canonicalize(colored1) == canonicalize(colored2)

        family_changed = [
            (c(Out()), bar(q)(Bulk(1))),
            (c(Bulk(1)), bar(q)(Bulk(2))),
            (c(Bulk(2)), bar(q)(In())),
        ]
        @test canonicalize(colored1) != canonicalize(family_changed)

        propagator_changed = [
            (c(Out()), bar(q)(Bulk(1))),
            (χc(Bulk(1)), bar(χc)(Bulk(2))),
            (c(Bulk(2)), bar(q)(In())),
        ]
        @test canonicalize(colored1) != canonicalize(propagator_changed)

        χ_up = FieldFamily{Boson}(:χ, FieldIndices(FieldIndex(:spin, 1)))
        χ_down = FieldFamily{Boson}(:χ, FieldIndices(FieldIndex(:spin, 2)))
        up_c, up_q = χ_up[Classical], χ_up[Quantum]
        down_c, down_q = χ_down[Classical], χ_down[Quantum]
        spin_up = [
            (c(Out()), bar(q)(Bulk(1))),
            (up_c(Bulk(1)), bar(up_q)(Bulk(2))),
            (c(Bulk(2)), bar(q)(In())),
        ]
        spin_down = [
            (c(Out()), bar(q)(Bulk(1))),
            (down_c(Bulk(1)), bar(down_q)(Bulk(2))),
            (c(Bulk(2)), bar(q)(In())),
        ]
        @test canonicalize(spin_up) != canonicalize(spin_down)
    end

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
        @inferred canonicalize(vs1)
    end
end
