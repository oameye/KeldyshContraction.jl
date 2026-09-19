using KeldyshContraction
using Combinatorics
import KeldyshContraction as KC

struct DirectSearchStats
    completed::Int
    raw_keys::Int
    cancelled_keys::Int
    filter_survivors::Int
end

struct GCSearchStats
    quotient_states::Int
    completed_states::Int
    transitions::Int
    merged_transitions::Int
    internal_canonicalization::Int
    pruned_transitions::Int
    filter_survivors::Int
end

function semantic_weights(pairings)
    weights = Dict{Any,Int}()
    for (pairing, topology, multiplicity) in pairings
        key = (pairing.contractions, topology)
        weight = Int(pairing.sign) * multiplicity
        weights[key] = get(weights, key, 0) + weight
    end
    filter!(pair -> !iszero(last(pair)), weights)
    return weights
end

function workload_terms(in_out, L, ::Val{O}) where {O}
    result = Vector{typeof(copy(in_out.args_nc))}()
    for coefficients in Combinatorics.multiexponents(length(L.lagrangian), O)
        idxs = KC.indices_from_counts(coefficients)
        mult = Combinatorics.multinomial(coefficients...)
        qmul = mult * prod(L(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))
        push!(result, copy((in_out * qmul).args_nc))
    end
    return result
end

function interaction_terms(::Type{S}, L, ::Val{O}) where {S<:KC.Statistics,O}
    products = KC.propagator_external_products(S, KC.propagator_fields(L, nothing)...)
    result = Vector{typeof(copy(first(products).args_nc))}()
    for in_out in products
        append!(result, workload_terms(in_out, L, Val(O)))
    end
    return result
end

function direct_search_stats(
    args_nc::Vector{KC.Field{S}}, ::Val{E}; regularise=true, _set_reg_to_zero=true
) where {S<:KC.Statistics,E}
    destroys, creates = KC.prepare_args(args_nc, Val(E))
    ps = map(KC.position, args_nc)
    candidates = KC.wick_candidates(
        destroys,
        creates,
        Val(E);
        regularise,
        _set_reg_to_zero,
        skip_external_pair=KC.has_in(ps) && KC.has_out(ps),
    )

    completed = 0
    matching_weights = Dict{NTuple{E,UInt8},Int}()
    KC.foreach_wick_matching(candidates, Val(E)) do contractions, permutation
        completed += 1
        key = KC.wick_matching_key(contractions, candidates, Val(E))
        weight = Int(KC.pairing_sign(S, permutation))
        matching_weights[key] = get(matching_weights, key, 0) + weight
        return nothing
    end

    cancelled = count(iszero, values(matching_weights))
    survivors = 0
    for (key, weight) in matching_weights
        iszero(weight) && continue
        contractions = KC.contractions_from_matching_key(key, candidates, KC.Contraction{S})
        KC.passes_wick_filters(contractions) && (survivors += 1)
    end
    return DirectSearchStats(completed, length(matching_weights), cancelled, survivors)
end

function gc_search_stats(
    args_nc::Vector{KC.Field{S}}, ::Val{E}; regularise=true, _set_reg_to_zero=true
) where {S<:KC.Statistics,E}
    problem, lookup, transport = KC._gc_build_wick_problem(
        args_nc, Val(E); regularise, _set_reg_to_zero
    )
    policy = KC._gc_causal_policy(lookup)
    completions, stats = KC.GC.generate_weighted_with_stats(problem; policy, transport)
    survivors = 0
    for completion in completions
        contractions = KC._gc_completion_contractions(completion, lookup, Val(E))
        KC.is_connected(contractions) && (survivors += 1)
    end
    return GCSearchStats(
        sum(stats.layer_states),
        isempty(stats.layer_states) ? 0 : last(stats.layer_states),
        stats.transitions,
        stats.merged_transitions,
        stats.canonicalization_calls,
        stats.pruned_transitions,
        survivors,
    )
end

function Base.:+(a::DirectSearchStats, b::DirectSearchStats)
    return DirectSearchStats(
        a.completed + b.completed,
        a.raw_keys + b.raw_keys,
        a.cancelled_keys + b.cancelled_keys,
        a.filter_survivors + b.filter_survivors,
    )
end

function Base.:+(a::GCSearchStats, b::GCSearchStats)
    return GCSearchStats(
        a.quotient_states + b.quotient_states,
        a.completed_states + b.completed_states,
        a.transitions + b.transitions,
        a.merged_transitions + b.merged_transitions,
        a.internal_canonicalization + b.internal_canonicalization,
        a.pruned_transitions + b.pruned_transitions,
        a.filter_survivors + b.filter_survivors,
    )
end

function direct_run(
    terms, ::Val{E}, ::Val{E2}; regularise=true, _set_reg_to_zero=true, simplify=true
) where {E,E2}
    output = 0
    for args_nc in terms
        output += length(
            KC._wick_contraction(
                args_nc, Val(E), Val(E2); regularise, _set_reg_to_zero, simplify
            ),
        )
    end
    return output
end

function gc_run(
    terms, ::Val{E}, ::Val{E2}; regularise=true, _set_reg_to_zero=true, simplify=true
) where {E,E2}
    output = 0
    for args_nc in terms
        output += length(
            KC._gc_wick_contraction(
                args_nc, Val(E), Val(E2); regularise, _set_reg_to_zero, simplify
            ),
        )
    end
    return output
end

function gc_build_all(terms, ::Val{E}; regularise=true, _set_reg_to_zero=true) where {E}
    return [
        KC._gc_build_wick_problem(args_nc, Val(E); regularise, _set_reg_to_zero) for
        args_nc in terms
    ]
end

function gc_generate_all(prepared)
    return [
        KC.GC.generate_weighted(
            first(x); policy=KC._gc_causal_policy(x[2]), transport=x[3]
        ) for x in prepared
    ]
end

function best_measurement(f, samples::Int)
    f()
    best_time = Inf
    best_bytes = typemax(Int)
    for _ in 1:samples
        measurement = @timed f()
        best_time = min(best_time, measurement.time)
        best_bytes = min(best_bytes, measurement.bytes)
    end
    return best_time, best_bytes
end

function certify_terms(
    terms, ::Val{E}, ::Val{E2}; regularise=true, _set_reg_to_zero=true, simplify=true
) where {E,E2}
    for args_nc in terms
        direct = KC._wick_contraction(
            args_nc, Val(E), Val(E2); regularise, _set_reg_to_zero, simplify
        )
        gc = KC._gc_wick_contraction(
            args_nc, Val(E), Val(E2); regularise, _set_reg_to_zero, simplify
        )
        semantic_weights(gc) == semantic_weights(direct) ||
            error("GC/direct semantic mismatch in production crossover benchmark")
    end
    return nothing
end

function benchmark_workload(name, ::Type{S}, L, ::Val{O}, ::Val{E}; samples=3) where {S,O,E}
    terms = interaction_terms(S, L, Val(O))
    E2 = KC.max_edges(O)
    regularise = KC.should_regularise(L.lagrangian)
    kwargs = (; regularise, _set_reg_to_zero=true, simplify=true)

    certify_terms(terms, Val(E), Val(E2); kwargs...)

    direct_stats = DirectSearchStats(0, 0, 0, 0)
    gc_stats = GCSearchStats(0, 0, 0, 0, 0, 0, 0)
    for args_nc in terms
        direct_stats += direct_search_stats(
            args_nc, Val(E); regularise, _set_reg_to_zero=true
        )
        gc_stats += gc_search_stats(args_nc, Val(E); regularise, _set_reg_to_zero=true)
    end

    direct = () -> direct_run(terms, Val(E), Val(E2); kwargs...)
    gc = () -> gc_run(terms, Val(E), Val(E2); kwargs...)
    direct_outputs = direct()
    gc_outputs = gc()
    direct_outputs == gc_outputs || error("GC/direct output-count mismatch")

    direct_time, direct_bytes = best_measurement(direct, samples)
    gc_time, gc_bytes = best_measurement(gc, samples)

    build = () -> gc_build_all(terms, Val(E); regularise, _set_reg_to_zero=true)
    build_time, build_bytes = best_measurement(build, samples)
    prepared = build()
    generate = () -> gc_generate_all(prepared)
    generate_time, generate_bytes = best_measurement(generate, samples)

    println(
        join(
            (
                name,
                string(S),
                O,
                E,
                length(terms),
                direct_stats.completed,
                direct_stats.raw_keys,
                direct_stats.cancelled_keys,
                direct_stats.filter_survivors,
                gc_stats.quotient_states,
                gc_stats.completed_states,
                gc_stats.transitions,
                gc_stats.merged_transitions,
                gc_stats.internal_canonicalization,
                gc_stats.pruned_transitions,
                gc_stats.filter_survivors,
                direct_outputs,
                direct_time,
                gc_time,
                gc_time / direct_time,
                direct_bytes,
                gc_bytes,
                gc_bytes / direct_bytes,
                build_time,
                build_bytes,
                generate_time,
                generate_bytes,
            ),
            '\t',
        ),
    )
    return nothing
end

@qfields bench_ϕ::Boson
c, q = bench_ϕ[Classical], bench_ϕ[Quantum]
boson_vertex = -(
    0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2)
)
L_b = InteractionLagrangian(boson_vertex, :g)

@qfields bench_ψ::Fermion
ψ₁, ψ₂ = bench_ψ[One], bench_ψ[Two]
fermion_vertex = ψ₁ * ψ₂ * bar(ψ₁) * bar(ψ₂)
L_f = InteractionLagrangian(fermion_vertex, :u)

println(
    "workload\tstatistics\torder\tedges\tterms\tdirect_completed\tdirect_raw_keys\t",
    "direct_cancelled_keys\tdirect_physical_canon\tgc_quotient_states\tgc_completed_states\t",
    "gc_transitions\tgc_merged_transitions\tgc_internal_canon\tgc_pruned_transitions\t",
    "gc_physical_canon\toutputs\tdirect_s\tgc_s\ttime_ratio\tdirect_bytes\tgc_bytes\t",
    "alloc_ratio\tgc_build_s\tgc_build_bytes\tgc_generate_s\tgc_generate_bytes",
)

benchmark_workload("boson_g2", Boson, L_b, Val(2), Val(5); samples=3)
benchmark_workload("boson_g3", Boson, L_b, Val(3), Val(7); samples=3)
benchmark_workload("boson_g4", Boson, L_b, Val(4), Val(9); samples=3)
benchmark_workload("fermion_u2", Fermion, L_f, Val(2), Val(5); samples=3)
benchmark_workload("fermion_u3", Fermion, L_f, Val(3), Val(7); samples=3)
