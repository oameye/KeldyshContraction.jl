using Test
using Combinatorics
using KeldyshContraction

import KeldyshContraction as KC

mutable struct ProductionCanonicalizationStats
    calls::Int
    oracle_equal::Int
    topology_equal::Int
    term_weight_checks::Int
end
ProductionCanonicalizationStats() = ProductionCanonicalizationStats(0, 0, 0, 0)

function add!(a::ProductionCanonicalizationStats, b::ProductionCanonicalizationStats)
    a.calls += b.calls
    a.oracle_equal += b.oracle_equal
    a.topology_equal += b.topology_equal
    a.term_weight_checks += b.term_weight_checks
    return a
end

function nauty_physical_canonicalize(vs::Vector{T}) where {T}
    isempty(vs) && return copy(vs)
    graph_positions = KC.canonicalization_positions(vs)
    physical_permutation, _, _, _ = KC.canonicalization_permutations(vs, graph_positions)
    mapping = KC.make_permutation_dict(physical_permutation, graph_positions, vs)
    return T[KC.relabel_bulk_positions(item, mapping) for item in vs]
end

function compare_fields(
    args_nc::Vector{KC.Field{S}},
    ::Val{E},
    ::Val{E2},
    scratch::KC.PhysicalCanonicalizationWorkspace;
    regularise=true,
    _set_reg_to_zero=false,
    simplify=false,
    collect_corpus=false,
) where {S<:KC.Statistics,E,E2}
    destroys, creates = KC.prepare_args(args_nc, Val(E))
    ps = map(KC.position, args_nc)
    skip = KC.has_in(ps) && KC.has_out(ps)
    candidates = KC.wick_candidates(
        destroys, creates, Val(E); regularise, _set_reg_to_zero, skip_external_pair=skip
    )

    matching_weights = Dict{NTuple{E,UInt8},Int}()
    KC.foreach_wick_matching(candidates, Val(E)) do contractions, permutation
        key = KC.wick_matching_key(contractions, candidates, Val(E))
        weight = Int(KC.pairing_sign(S, permutation))
        matching_weights[key] = get(matching_weights, key, 0) + weight
        return nothing
    end

    nauty_weights = Dict{Any,Int}()
    production_weights = Dict{Any,Int}()
    corpus = Vector{Vector{KC.Contraction{S}}}()
    stats = ProductionCanonicalizationStats()

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

        nauty = nauty_physical_canonicalize(final_contractions)
        production = KC.canonicalize(final_contractions, scratch)

        nauty_topology = KC.legacy_topology(nauty, Val(E2))
        production_topology = KC.legacy_topology(production, Val(E2))
        topology_equal = production_topology == nauty_topology
        stats.topology_equal += Int(topology_equal)
        @test topology_equal

        # Normalize the production representative with the independent old Nauty path only
        # for certification. GC and Nauty intentionally use different private canonical ranks.
        production_oracle = nauty_physical_canonicalize(production)
        nauty_key = KC.sorted_wick_key(nauty, Val(E))
        production_key = KC.sorted_wick_key(production_oracle, Val(E))
        oracle_equal = production_key == nauty_key
        stats.oracle_equal += Int(oracle_equal)
        @test oracle_equal

        nauty_weights[nauty_key] = get(nauty_weights, nauty_key, 0) + final_weight
        production_weights[production_key] =
            get(production_weights, production_key, 0) + final_weight
    end

    filter!((kv) -> !iszero(kv[2]), nauty_weights)
    filter!((kv) -> !iszero(kv[2]), production_weights)
    @test production_weights == nauty_weights
    stats.term_weight_checks += 1
    return stats, corpus
end

function compare_component(
    in_out::KC.QMul{CI,S},
    L::KC.InteractionLagrangian{CL,S},
    ::Val{O},
    ::Val{E},
    scratch::KC.PhysicalCanonicalizationWorkspace;
    simplify=true,
    _set_reg_to_zero=true,
    collect_corpus=false,
) where {CI<:Number,CL<:Number,S<:KC.Statistics,O,E}
    total = ProductionCanonicalizationStats()
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

statistics_type(::KC.InteractionLagrangian{C,S}) where {C,S} = S

function compare_interaction(
    label,
    L::KC.InteractionLagrangian,
    order::Int,
    edges::Int,
    scratch::KC.PhysicalCanonicalizationWorkspace;
    simplify=true,
    _set_reg_to_zero=true,
    collect_corpus=false,
)
    S = statistics_type(L)
    fields = KC.propagator_fields(L, nothing)
    products = KC.propagator_external_products(S, fields...)
    total = ProductionCanonicalizationStats()
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
            "PRODUCTION_SEMANTICS\t",
            label,
            "\t",
            component,
            "\tcalls=",
            stats.calls,
            "\toracle_equal=",
            stats.oracle_equal,
            "\ttopology_equal=",
            stats.topology_equal,
            "\tterm_weights=",
            stats.term_weight_checks,
        )
    end

    @test total.oracle_equal == total.calls
    @test total.topology_equal == total.calls
    return total, corpus
end

function best_repeated_time(f, repetitions; samples=7)
    best = Inf
    for _ in 1:samples
        elapsed = @elapsed for _ in 1:repetitions
            f()
        end
        best = min(best, elapsed)
    end
    return best / repetitions
end

function benchmark_corpus(label, corpus, scratch)
    isempty(corpus) && return nothing
    production_batch = () -> begin
        for contractions in corpus
            KC.canonicalize(contractions, scratch)
        end
        nothing
    end
    nauty_batch = () -> begin
        for contractions in corpus
            nauty_physical_canonicalize(contractions)
        end
        nothing
    end

    production_batch()
    nauty_batch()
    repetitions = max(1, cld(50_000, length(corpus)))
    production_time = best_repeated_time(production_batch, repetitions)
    nauty_time = best_repeated_time(nauty_batch, repetitions)
    production_bytes = @allocated production_batch()
    nauty_bytes = @allocated nauty_batch()
    time_ratio = production_time / nauty_time
    alloc_ratio = production_bytes / nauty_bytes

    println(
        "PRODUCTION_CORPUS\t",
        label,
        "\tcalls=",
        length(corpus),
        "\ttiming_repetitions=",
        repetitions,
        "\ttimed_calls=",
        repetitions * length(corpus),
        "\ttime_ratio=",
        round(time_ratio; digits=3),
        "\talloc_ratio=",
        round(alloc_ratio; digits=3),
        "\tgc_ms=",
        round(production_time * 1e3; digits=3),
        "\tnauty_ms=",
        round(nauty_time * 1e3; digits=3),
        "\tgc_bytes=",
        production_bytes,
        "\tnauty_bytes=",
        nauty_bytes,
    )
    @test time_ratio <= 1.0
    @test alloc_ratio <= 1.0
    return nothing
end

scratch = KC.PhysicalCanonicalizationWorkspace(32)

@qfields ϕwick::Boson
c, q = ϕwick[Classical], ϕwick[Quantum]

elastic = -(0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2))
L_g = InteractionLagrangian(elastic, :g)

loss =
    0.5 *
    bar(c) *
    bar(q) *
    (
        c(KC.Regularisation.Minus) * c(KC.Regularisation.Minus) +
        q(KC.Regularisation.Minus) * q(KC.Regularisation.Minus)
    ) -
    0.5 *
    c(KC.Regularisation.Plus) *
    q(KC.Regularisation.Plus) *
    (bar(c) * bar(c) + bar(q) * bar(q)) +
    bar(c) *
    bar(q) *
    (
        c(KC.Regularisation.Plus) * q(KC.Regularisation.Plus) +
        c(KC.Regularisation.Minus) * q(KC.Regularisation.Minus)
    )
L_γ = InteractionLagrangian(loss, :γ)

@testset "production GC physical canonicalization on real Wick workloads" begin
    total_calls = 0
    stats, g2_corpus = compare_interaction(
        "boson_g2", L_g, 2, 5, scratch; collect_corpus=true
    )
    total_calls += stats.calls
    stats, gamma2_corpus = compare_interaction(
        "boson_gamma2",
        L_γ,
        2,
        5,
        scratch;
        simplify=true,
        _set_reg_to_zero=true,
        collect_corpus=true,
    )
    total_calls += stats.calls
    stats, g3_corpus = compare_interaction(
        "boson_g3", L_g, 3, 7, scratch; collect_corpus=true
    )
    total_calls += stats.calls

    @qfields ψwick::Fermion
    ψ₁, ψ₂ = ψwick[One], ψwick[Two]
    fermion_vertex = ψ₁ * ψ₂ * bar(ψ₁) * bar(ψ₂)
    L_f = InteractionLagrangian(fermion_vertex, :u)
    stats, f2_corpus = compare_interaction(
        "fermion_quartic2", L_f, 2, 5, scratch; simplify=false, collect_corpus=true
    )
    total_calls += stats.calls
    stats, f3_corpus = compare_interaction(
        "fermion_quartic3", L_f, 3, 7, scratch; simplify=false, collect_corpus=true
    )
    total_calls += stats.calls

    ∂xψ₂ = partial(ψ₂, :x)
    derivative_vertex = ψ₁ * ∂xψ₂ * bar(ψ₁) * bar(∂xψ₂)
    L_p = InteractionLagrangian(derivative_vertex, :γ)
    stats, fp2_corpus = compare_interaction(
        "fermion_derivative2", L_p, 2, 5, scratch; simplify=false, collect_corpus=true
    )
    total_calls += stats.calls

    @test total_calls == 14_347

    println()
    println("production GC vs independent old-Nauty real-corpus benchmark")
    benchmark_corpus("boson_g2", g2_corpus, scratch)
    benchmark_corpus("boson_gamma2", gamma2_corpus, scratch)
    benchmark_corpus("boson_g3", g3_corpus, scratch)
    benchmark_corpus("fermion_quartic2", f2_corpus, scratch)
    benchmark_corpus("fermion_quartic3", f3_corpus, scratch)
    benchmark_corpus("fermion_derivative2", fp2_corpus, scratch)
end
