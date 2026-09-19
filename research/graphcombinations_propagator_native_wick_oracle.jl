using Test
using Combinatorics
using KeldyshContraction

import KeldyshContraction as KC

include(joinpath(@__DIR__, "graphcombinations_propagator_gc_adapter.jl"))
include(joinpath(@__DIR__, "graphcombinations_propagator_native_adapter.jl"))

mutable struct NativeWickCanonicalizationStats
    calls::Int
    hybrid_oracle_equal::Int
    native_oracle_equal::Int
    hybrid_topology_equal::Int
    native_topology_equal::Int
    term_weight_checks::Int
end
NativeWickCanonicalizationStats() = NativeWickCanonicalizationStats(0, 0, 0, 0, 0, 0)

function add!(a::NativeWickCanonicalizationStats, b::NativeWickCanonicalizationStats)
    a.calls += b.calls
    a.hybrid_oracle_equal += b.hybrid_oracle_equal
    a.native_oracle_equal += b.native_oracle_equal
    a.hybrid_topology_equal += b.hybrid_topology_equal
    a.native_topology_equal += b.native_topology_equal
    a.term_weight_checks += b.term_weight_checks
    return a
end

function compare_native_fields(
    args_nc::Vector{KC.Field{S}},
    ::Val{E},
    ::Val{E2},
    scratch::GCNativePhysicalCanonicalizationWorkspace;
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
    hybrid_weights = Dict{Any,Int}()
    native_weights = Dict{Any,Int}()
    corpus = Vector{Vector{KC.Contraction{S}}}()
    stats = NativeWickCanonicalizationStats()

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
        hybrid = gc_native_hybrid_physical_canonicalize!(scratch, final_contractions)
        native = gc_native_physical_canonicalize!(scratch, final_contractions)

        nauty_topology = KC.legacy_topology(nauty, Val(E2))
        hybrid_topology = KC.legacy_topology(hybrid, Val(E2))
        native_topology = KC.legacy_topology(native, Val(E2))
        hybrid_topology_equal = hybrid_topology == nauty_topology
        native_topology_equal = native_topology == nauty_topology
        stats.hybrid_topology_equal += Int(hybrid_topology_equal)
        stats.native_topology_equal += Int(native_topology_equal)
        @test hybrid_topology_equal
        @test native_topology_equal

        nauty_key = KC.sorted_wick_key(nauty, Val(E))
        hybrid_key = KC.sorted_wick_key(KC.canonicalize(hybrid), Val(E))
        native_key = KC.sorted_wick_key(KC.canonicalize(native), Val(E))
        hybrid_equal = hybrid_key == nauty_key
        native_equal = native_key == nauty_key
        stats.hybrid_oracle_equal += Int(hybrid_equal)
        stats.native_oracle_equal += Int(native_equal)
        @test hybrid_equal
        @test native_equal

        baseline_weights[nauty_key] = get(baseline_weights, nauty_key, 0) + final_weight
        hybrid_weights[hybrid_key] = get(hybrid_weights, hybrid_key, 0) + final_weight
        native_weights[native_key] = get(native_weights, native_key, 0) + final_weight
    end

    filter!((kv) -> !iszero(kv[2]), baseline_weights)
    filter!((kv) -> !iszero(kv[2]), hybrid_weights)
    filter!((kv) -> !iszero(kv[2]), native_weights)
    @test hybrid_weights == baseline_weights
    @test native_weights == baseline_weights
    stats.term_weight_checks += 1
    return stats, corpus
end

function compare_native_component(
    in_out::KC.QMul{CI,S},
    L::KC.InteractionLagrangian{CL,S},
    ::Val{O},
    ::Val{E},
    scratch::GCNativePhysicalCanonicalizationWorkspace;
    simplify=true,
    _set_reg_to_zero=true,
    collect_corpus=false,
) where {CI<:Number,CL<:Number,S<:KC.Statistics,O,E}
    total = NativeWickCanonicalizationStats()
    corpus = Vector{Vector{KC.Contraction{S}}}()
    regularise = KC.should_regularise(L.lagrangian)
    l = length(L.lagrangian)

    for coefficients in Combinatorics.multiexponents(l, O)
        idxs = KC.indices_from_counts(coefficients)
        mult = Combinatorics.multinomial(coefficients...)
        qmul = mult * prod(L(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))
        term = in_out * qmul
        stats, term_corpus = compare_native_fields(
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

function native_statistics_type(::KC.InteractionLagrangian{C,S}) where {C,S}
    return S
end

function compare_native_interaction(
    label,
    L::KC.InteractionLagrangian,
    order::Int,
    edges::Int,
    scratch::GCNativePhysicalCanonicalizationWorkspace;
    simplify=true,
    _set_reg_to_zero=true,
    collect_corpus=false,
)
    S = native_statistics_type(L)
    fields = KC.propagator_fields(L, nothing)
    products = KC.propagator_external_products(S, fields...)
    total = NativeWickCanonicalizationStats()
    corpus = Vector{Vector{KC.Contraction{S}}}()

    for (component, in_out) in zip(("K", "R", "A"), products)
        stats, component_corpus = compare_native_component(
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
            "NATIVE_SEMANTICS\t", label, "\t", component,
            "\tcalls=", stats.calls,
            "\thybrid_oracle=", stats.hybrid_oracle_equal,
            "\tnative_oracle=", stats.native_oracle_equal,
            "\thybrid_topology=", stats.hybrid_topology_equal,
            "\tnative_topology=", stats.native_topology_equal,
            "\tterm_weights=", stats.term_weight_checks,
        )
    end

    @test total.hybrid_oracle_equal == total.calls
    @test total.native_oracle_equal == total.calls
    @test total.hybrid_topology_equal == total.calls
    @test total.native_topology_equal == total.calls
    return total, corpus
end

function native_best_corpus_time(f; samples=5)
    best = Inf
    for _ in 1:samples
        elapsed = @elapsed f()
        best = min(best, elapsed)
    end
    return best
end

function benchmark_native_corpus(label, corpus, gadget_scratch, native_scratch)
    isempty(corpus) && return nothing
    gadget_batch = () -> begin
        for contractions in corpus
            gc_physical_canonicalize!(gadget_scratch, contractions)
        end
        nothing
    end
    native_batch = () -> begin
        for contractions in corpus
            gc_native_physical_canonicalize!(native_scratch, contractions)
        end
        nothing
    end
    hybrid_batch = () -> begin
        for contractions in corpus
            gc_native_hybrid_physical_canonicalize!(native_scratch, contractions)
        end
        nothing
    end
    nauty_batch = () -> begin
        for contractions in corpus
            KC.canonicalize(contractions)
        end
        nothing
    end

    gadget_batch()
    native_batch()
    hybrid_batch()
    nauty_batch()
    gadget_time = native_best_corpus_time(gadget_batch)
    native_time = native_best_corpus_time(native_batch)
    hybrid_time = native_best_corpus_time(hybrid_batch)
    nauty_time = native_best_corpus_time(nauty_batch)
    gadget_bytes = @allocated gadget_batch()
    native_bytes = @allocated native_batch()
    hybrid_bytes = @allocated hybrid_batch()
    nauty_bytes = @allocated nauty_batch()
    println(
        "NATIVE_CORPUS\t", label,
        "\tcalls=", length(corpus),
        "\tgadget_over_nauty=", round(gadget_time / nauty_time; digits=3),
        "\tnative_over_nauty=", round(native_time / nauty_time; digits=3),
        "\thybrid_over_nauty=", round(hybrid_time / nauty_time; digits=3),
        "\tnative_over_gadget=", round(native_time / gadget_time; digits=3),
        "\thybrid_over_gadget=", round(hybrid_time / gadget_time; digits=3),
        "\tgadget_ms=", round(gadget_time * 1e3; digits=3),
        "\tnative_ms=", round(native_time * 1e3; digits=3),
        "\thybrid_ms=", round(hybrid_time * 1e3; digits=3),
        "\tnauty_ms=", round(nauty_time * 1e3; digits=3),
        "\tgadget_bytes=", gadget_bytes,
        "\tnative_bytes=", native_bytes,
        "\thybrid_bytes=", hybrid_bytes,
        "\tnauty_bytes=", nauty_bytes,
    )
    return nothing
end

gadget_scratch = GCPhysicalCanonicalizationWorkspace(32)
native_scratch = GCNativePhysicalCanonicalizationWorkspace(32)

@qfields ϕnative::Boson
c, q = ϕnative[Classical], ϕnative[Quantum]

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

@testset "native GC physical canonicalization on real Wick workloads" begin
    totals = NativeWickCanonicalizationStats()
    stats, g2_corpus = compare_native_interaction(
        "boson_g2", L_g, 2, 5, native_scratch; collect_corpus=true
    )
    add!(totals, stats)
    stats, gamma2_corpus = compare_native_interaction(
        "boson_gamma2", L_γ, 2, 5, native_scratch;
        simplify=true, _set_reg_to_zero=true, collect_corpus=true,
    )
    add!(totals, stats)
    stats, g3_corpus = compare_native_interaction(
        "boson_g3", L_g, 3, 7, native_scratch; collect_corpus=true
    )
    add!(totals, stats)

    @qfields ψnative::Fermion
    ψ₁, ψ₂ = ψnative[One], ψnative[Two]
    fermion_vertex = ψ₁ * ψ₂ * bar(ψ₁) * bar(ψ₂)
    L_f = InteractionLagrangian(fermion_vertex, :u)
    stats, f2_corpus = compare_native_interaction(
        "fermion_quartic2", L_f, 2, 5, native_scratch;
        simplify=false, collect_corpus=true,
    )
    add!(totals, stats)
    stats, f3_corpus = compare_native_interaction(
        "fermion_quartic3", L_f, 3, 7, native_scratch;
        simplify=false, collect_corpus=true,
    )
    add!(totals, stats)

    ∂xψ₂ = partial(ψ₂, :x)
    derivative_vertex = ψ₁ * ∂xψ₂ * bar(ψ₁) * bar(∂xψ₂)
    L_p = InteractionLagrangian(derivative_vertex, :γ)
    stats, fp2_corpus = compare_native_interaction(
        "fermion_derivative2", L_p, 2, 5, native_scratch;
        simplify=false, collect_corpus=true,
    )
    add!(totals, stats)

    @test totals.calls == 28_868
    @test totals.hybrid_oracle_equal == totals.calls
    @test totals.native_oracle_equal == totals.calls
    @test totals.hybrid_topology_equal == totals.calls
    @test totals.native_topology_equal == totals.calls
    println(
        "NATIVE_TOTAL\tcalls=", totals.calls,
        "\thybrid_oracle=", totals.hybrid_oracle_equal,
        "\tnative_oracle=", totals.native_oracle_equal,
        "\thybrid_topology=", totals.hybrid_topology_equal,
        "\tnative_topology=", totals.native_topology_equal,
    )

    println()
    println("native relation real-corpus canonicalization benchmark")
    benchmark_native_corpus("boson_g2", g2_corpus, gadget_scratch, native_scratch)
    benchmark_native_corpus("boson_gamma2", gamma2_corpus, gadget_scratch, native_scratch)
    benchmark_native_corpus("boson_g3", g3_corpus, gadget_scratch, native_scratch)
    benchmark_native_corpus("fermion_quartic2", f2_corpus, gadget_scratch, native_scratch)
    benchmark_native_corpus("fermion_quartic3", f3_corpus, gadget_scratch, native_scratch)
    benchmark_native_corpus("fermion_derivative2", fp2_corpus, gadget_scratch, native_scratch)
end
