using Test
using Combinatorics
using KeldyshContraction

import KeldyshContraction as KC

include(joinpath(@__DIR__, "graphcombinations_propagator_gc_adapter.jl"))

mutable struct WickCanonicalizationStats
    calls::Int
    direct_equal::Int
    oracle_equal::Int
    topology_equal::Int
    term_weight_checks::Int
end
WickCanonicalizationStats() = WickCanonicalizationStats(0, 0, 0, 0, 0)

function add!(a::WickCanonicalizationStats, b::WickCanonicalizationStats)
    a.calls += b.calls
    a.direct_equal += b.direct_equal
    a.oracle_equal += b.oracle_equal
    a.topology_equal += b.topology_equal
    a.term_weight_checks += b.term_weight_checks
    return a
end

function compare_fields(
    args_nc::Vector{KC.Field{S}},
    ::Val{E},
    ::Val{E2},
    scratch::GCPhysicalCanonicalizationWorkspace;
    regularise=true,
    _set_reg_to_zero=false,
    simplify=false,
    collect_corpus=false,
) where {S<:KC.Statistics,E,E2}
    destroys, creates = KC.prepare_args(args_nc, Val(E))
    ps = map(KC.position, args_nc)
    skip = KC.has_in(ps) && KC.has_out(ps)
    candidates = KC.wick_candidates(
        destroys,
        creates,
        Val(E);
        regularise,
        _set_reg_to_zero,
        skip_external_pair=skip,
    )

    matching_weights = Dict{NTuple{E,UInt8},Int}()
    KC.foreach_wick_matching(candidates, Val(E)) do contractions, permutation
        key = KC.wick_matching_key(contractions, candidates, Val(E))
        weight = Int(KC.pairing_sign(S, permutation))
        matching_weights[key] = get(matching_weights, key, 0) + weight
        return nothing
    end

    baseline_weights = Dict{Any,Int}()
    gc_oracle_weights = Dict{Any,Int}()
    corpus = Vector{Vector{KC.Contraction{S}}}()
    stats = WickCanonicalizationStats()

    for (key, weight) in matching_weights
        iszero(weight) && continue
        contractions = KC.contractions_from_matching_key(key, candidates, KC.Contraction{S})
        KC.passes_wick_filters(contractions) || continue

        final_contractions, simplification_sign = if simplify
            KC.advanced_to_retarded(contractions, 1)
        else
            contractions, 1
        end
        final_weight = weight * Int(simplification_sign)
        stats.calls += 1
        collect_corpus && push!(corpus, copy(final_contractions))

        nauty = KC.canonicalize(final_contractions)
        gc = gc_physical_canonicalize!(scratch, final_contractions)

        stats.direct_equal += Int(gc == nauty)

        nauty_topology = KC.legacy_topology(nauty, Val(E2))
        gc_topology = KC.legacy_topology(gc, Val(E2))
        topology_equal = gc_topology == nauty_topology
        stats.topology_equal += Int(topology_equal)
        @test topology_equal

        # Normalize the GC representative with the independent current Nauty path only for
        # certification. This checks physical equivalence without requiring identical private
        # canonical-label conventions.
        gc_oracle = KC.canonicalize(gc)
        nauty_key = KC.sorted_wick_key(nauty, Val(E))
        gc_oracle_key = KC.sorted_wick_key(gc_oracle, Val(E))
        oracle_equal = gc_oracle_key == nauty_key
        stats.oracle_equal += Int(oracle_equal)
        @test oracle_equal

        baseline_weights[nauty_key] = get(baseline_weights, nauty_key, 0) + final_weight
        gc_oracle_weights[gc_oracle_key] =
            get(gc_oracle_weights, gc_oracle_key, 0) + final_weight
    end

    filter!((kv) -> !iszero(kv[2]), baseline_weights)
    filter!((kv) -> !iszero(kv[2]), gc_oracle_weights)
    @test gc_oracle_weights == baseline_weights
    stats.term_weight_checks += 1
    return stats, corpus
end

function compare_component(
    in_out::KC.QMul{CI,S},
    L::KC.InteractionLagrangian{CL,S},
    ::Val{O},
    ::Val{E},
    scratch::GCPhysicalCanonicalizationWorkspace;
    simplify=true,
    _set_reg_to_zero=true,
    collect_corpus=false,
) where {CI<:Number,CL<:Number,S<:KC.Statistics,O,E}
    total = WickCanonicalizationStats()
    corpus = Vector{Vector{KC.Contraction{S}}}()
    regularise = KC.should_regularise(L.lagrangian)
    l = length(L.lagrangian)

    for coefficients in Combinatorics.multiexponents(l, O)
        idxs = KC.indices_from_counts(coefficients)
        mult = Combinatorics.multinomial(coefficients...)
        qmul = mult * prod(L(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))
        term = in_out * qmul
        stats, term_corpus = compare_fields(
            term.args_nc,
            Val(E),
            Val(KC.max_edges(O)),
            scratch;
            regularise,
            _set_reg_to_zero,
            simplify,
            collect_corpus,
        )
        add!(total, stats)
        collect_corpus && append!(corpus, term_corpus)
    end
    return total, corpus
end

function statistics_type(::KC.InteractionLagrangian{C,S}) where {C,S}
    return S
end

function compare_interaction(
    label,
    L::KC.InteractionLagrangian,
    order::Int,
    edges::Int,
    scratch::GCPhysicalCanonicalizationWorkspace;
    simplify=true,
    _set_reg_to_zero=true,
    collect_corpus=false,
)
    S = statistics_type(L)
    fields = KC.propagator_fields(L, nothing)
    products = KC.propagator_external_products(S, fields...)
    total = WickCanonicalizationStats()
    corpus = Vector{Vector{KC.Contraction{S}}}()

    for (component, in_out) in zip(("K", "R", "A"), products)
        stats, component_corpus = compare_component(
            in_out,
            L,
            Val(order),
            Val(edges),
            scratch;
            simplify,
            _set_reg_to_zero,
            collect_corpus,
        )
        add!(total, stats)
        collect_corpus && append!(corpus, component_corpus)
        println(
            "SEMANTICS\t", label, "\t", component,
            "\tcalls=", stats.calls,
            "\tdirect_equal=", stats.direct_equal,
            "\toracle_equal=", stats.oracle_equal,
            "\ttopology_equal=", stats.topology_equal,
            "\tterm_weights=", stats.term_weight_checks,
        )
    end

    @test total.oracle_equal == total.calls
    @test total.topology_equal == total.calls
    println(
        "SEMANTICS\t", label, "\ttotal",
        "\tcalls=", total.calls,
        "\tdirect_equal=", total.direct_equal,
        "\toracle_equal=", total.oracle_equal,
        "\ttopology_equal=", total.topology_equal,
        "\tterm_weights=", total.term_weight_checks,
    )
    return total, corpus
end

function best_corpus_time(f; samples=5)
    best = Inf
    for _ in 1:samples
        elapsed = @elapsed f()
        best = min(best, elapsed)
    end
    return best
end

function benchmark_corpus(label, corpus, scratch)
    isempty(corpus) && return nothing
    gc_batch = () -> begin
        for contractions in corpus
            gc_physical_canonicalize!(scratch, contractions)
        end
        nothing
    end
    nauty_batch = () -> begin
        for contractions in corpus
            KC.canonicalize(contractions)
        end
        nothing
    end

    gc_batch()
    nauty_batch()
    gc_time = best_corpus_time(gc_batch)
    nauty_time = best_corpus_time(nauty_batch)
    gc_bytes = @allocated gc_batch()
    nauty_bytes = @allocated nauty_batch()
    println(
        "CORPUS\t", label,
        "\tcalls=", length(corpus),
        "\ttime_ratio=", round(gc_time / nauty_time; digits=3),
        "\talloc_ratio=", round(gc_bytes / nauty_bytes; digits=3),
        "\tgc_ms=", round(gc_time * 1e3; digits=3),
        "\tnauty_ms=", round(nauty_time * 1e3; digits=3),
        "\tgc_bytes=", gc_bytes,
        "\tnauty_bytes=", nauty_bytes,
    )
    return nothing
end

scratch = GCPhysicalCanonicalizationWorkspace(32)

@qfields ϕwick::Boson
c, q = ϕwick[Classical], ϕwick[Quantum]

elastic = -(
    0.5 * (c^2 + q^2) * bar(c) * bar(q) +
    0.5 * c * q * (bar(c)^2 + bar(q)^2)
)
L_g = InteractionLagrangian(elastic, :g)

loss =
    0.5 * bar(c) * bar(q) *
    (c(KC.Regularisation.Minus) * c(KC.Regularisation.Minus) +
     q(KC.Regularisation.Minus) * q(KC.Regularisation.Minus)) -
    0.5 * c(KC.Regularisation.Plus) * q(KC.Regularisation.Plus) *
    (bar(c) * bar(c) + bar(q) * bar(q)) +
    bar(c) * bar(q) *
    (c(KC.Regularisation.Plus) * q(KC.Regularisation.Plus) +
     c(KC.Regularisation.Minus) * q(KC.Regularisation.Minus))
L_γ = InteractionLagrangian(loss, :γ)

@testset "GC physical canonicalization on real Wick workloads" begin
    _, g2_corpus = compare_interaction("boson_g2", L_g, 2, 5, scratch; collect_corpus=true)
    _, gamma2_corpus = compare_interaction(
        "boson_gamma2", L_γ, 2, 5, scratch;
        simplify=true, _set_reg_to_zero=true, collect_corpus=true,
    )
    _, g3_corpus = compare_interaction("boson_g3", L_g, 3, 7, scratch; collect_corpus=true)

    @qfields ψwick::Fermion
    ψ₁, ψ₂ = ψwick[One], ψwick[Two]
    fermion_vertex = ψ₁ * ψ₂ * bar(ψ₁) * bar(ψ₂)
    L_f = InteractionLagrangian(fermion_vertex, :u)
    _, f2_corpus = compare_interaction(
        "fermion_quartic2", L_f, 2, 5, scratch; simplify=false, collect_corpus=true,
    )
    _, f3_corpus = compare_interaction(
        "fermion_quartic3", L_f, 3, 7, scratch; simplify=false, collect_corpus=true,
    )

    ∂xψ₂ = partial(ψ₂, :x)
    derivative_vertex = ψ₁ * ∂xψ₂ * bar(ψ₁) * bar(∂xψ₂)
    L_p = InteractionLagrangian(derivative_vertex, :γ)
    _, fp2_corpus = compare_interaction(
        "fermion_derivative2", L_p, 2, 5, scratch; simplify=false, collect_corpus=true,
    )

    println()
    println("production-shaped real-corpus canonicalization benchmark")
    benchmark_corpus("boson_g2", g2_corpus, scratch)
    benchmark_corpus("boson_gamma2", gamma2_corpus, scratch)
    benchmark_corpus("boson_g3", g3_corpus, scratch)
    benchmark_corpus("fermion_quartic2", f2_corpus, scratch)
    benchmark_corpus("fermion_quartic3", f3_corpus, scratch)
    benchmark_corpus("fermion_derivative2", fp2_corpus, scratch)
end
