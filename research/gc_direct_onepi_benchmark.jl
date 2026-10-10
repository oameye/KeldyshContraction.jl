using KeldyshContraction
using Combinatorics
import KeldyshContraction as KC

struct OnePIAggregateStats
    layer_states::Int
    completed_states::Int
    transitions::Int
    canonicalization_calls::Int
    pruned_transitions::Int
    outputs::Int
end

function Base.:+(a::OnePIAggregateStats, b::OnePIAggregateStats)
    return OnePIAggregateStats(
        a.layer_states + b.layer_states,
        a.completed_states + b.completed_states,
        a.transitions + b.transitions,
        a.canonicalization_calls + b.canonicalization_calls,
        a.pruned_transitions + b.pruned_transitions,
        a.outputs + b.outputs,
    )
end

function onepi_semantic_weights(pairings)
    weights = Dict{Any,Int}()
    for (pairing, topology, multiplicity) in pairings
        key = (pairing.contractions, topology)
        weight = Int(pairing.sign) * multiplicity
        weights[key] = get(weights, key, 0) + weight
    end
    filter!(pair -> !iszero(last(pair)), weights)
    return weights
end

function onepi_workload_terms(in_out, L, ::Val{O}) where {O}
    result = Vector{typeof(copy(in_out.args_nc))}()
    for coefficients in Combinatorics.multiexponents(length(L.lagrangian), O)
        idxs = KC.indices_from_counts(coefficients)
        mult = Combinatorics.multinomial(coefficients...)
        qmul = mult * prod(L(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))
        push!(result, copy((in_out * qmul).args_nc))
    end
    return result
end

function onepi_interaction_terms(::Type{S}, L, ::Val{O}) where {S<:KC.Statistics,O}
    products = KC.propagator_external_products(S, KC.propagator_fields(L, nothing)...)
    result = Vector{typeof(copy(first(products).args_nc))}()
    for in_out in products
        append!(result, onepi_workload_terms(in_out, L, Val(O)))
    end
    return result
end

function onepi_pairings(
    args_nc, ::Val{E}, ::Val{E2}; onepi_pruning, regularise, simplify
) where {E,E2}
    return KC._gc_wick_contraction_with_stats(
        args_nc,
        Val(E),
        Val(E2);
        regularise,
        _set_reg_to_zero=true,
        simplify,
        causal_pruning=true,
        onepi_pruning,
    )
end

function baseline_onepi(pairings)
    return [entry for entry in pairings if KC.is_irreducible(first(entry).contractions)]
end

function aggregate_onepi_stats(stats, outputs)
    return OnePIAggregateStats(
        sum(stats.layer_states),
        isempty(stats.layer_states) ? 0 : last(stats.layer_states),
        stats.transitions,
        stats.canonicalization_calls,
        stats.pruned_transitions,
        outputs,
    )
end

function run_onepi_corpus(
    terms, ::Val{E}, ::Val{E2}; onepi_pruning, regularise, simplify
) where {E,E2}
    aggregate = OnePIAggregateStats(0, 0, 0, 0, 0, 0)
    for args_nc in terms
        pairings, stats = onepi_pairings(
            args_nc, Val(E), Val(E2); onepi_pruning, regularise, simplify
        )
        outputs = onepi_pruning ? pairings : baseline_onepi(pairings)
        aggregate += aggregate_onepi_stats(stats, length(outputs))
    end
    return aggregate
end

function certify_onepi_corpus(terms, ::Val{E}, ::Val{E2}; regularise, simplify) where {E,E2}
    for args_nc in terms
        ordinary, _ = onepi_pairings(
            args_nc, Val(E), Val(E2); onepi_pruning=false, regularise, simplify
        )
        expected = baseline_onepi(ordinary)
        direct, _ = onepi_pairings(
            args_nc, Val(E), Val(E2); onepi_pruning=true, regularise, simplify
        )
        onepi_semantic_weights(direct) == onepi_semantic_weights(expected) ||
            error("direct 1PI semantic mismatch")
    end
    return nothing
end

function best_onepi_measurement(f; samples=2)
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

function benchmark_onepi_workload(
    name, ::Type{S}, L, ::Val{O}, ::Val{E}; simplify=true, samples=2
) where {S,O,E}
    terms = onepi_interaction_terms(S, L, Val(O))
    E2 = KC.max_edges(O)
    regularise = KC.should_regularise(L.lagrangian)

    certify_onepi_corpus(terms, Val(E), Val(E2); regularise, simplify)

    baseline =
        () -> run_onepi_corpus(
            terms, Val(E), Val(E2); onepi_pruning=false, regularise, simplify
        )
    direct =
        () -> run_onepi_corpus(
            terms, Val(E), Val(E2); onepi_pruning=true, regularise, simplify
        )

    baseline_stats = baseline()
    direct_stats = direct()
    baseline_stats.outputs == direct_stats.outputs ||
        error("direct 1PI output-count mismatch")

    baseline_time, baseline_bytes = best_onepi_measurement(baseline; samples)
    direct_time, direct_bytes = best_onepi_measurement(direct; samples)

    println(
        join(
            (
                name,
                string(S),
                O,
                E,
                length(terms),
                baseline_stats.layer_states,
                direct_stats.layer_states,
                baseline_stats.completed_states,
                direct_stats.completed_states,
                baseline_stats.transitions,
                direct_stats.transitions,
                baseline_stats.canonicalization_calls,
                direct_stats.canonicalization_calls,
                baseline_stats.pruned_transitions,
                direct_stats.pruned_transitions,
                direct_stats.pruned_transitions - baseline_stats.pruned_transitions,
                direct_stats.outputs,
                baseline_time,
                direct_time,
                direct_time / baseline_time,
                baseline_bytes,
                direct_bytes,
                direct_bytes / baseline_bytes,
            ),
            '\t',
        ),
    )
    return nothing
end

@qfields onepi_ϕ::Boson
c, q = onepi_ϕ[Classical], onepi_ϕ[Quantum]
boson_vertex = -(0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2))
L_b = InteractionLagrangian(boson_vertex, :g)

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

@qfields onepi_ψ::Fermion
ψ₁, ψ₂ = onepi_ψ[One], onepi_ψ[Two]
fermion_vertex = ψ₁ * ψ₂ * bar(ψ₁) * bar(ψ₂)
L_f = InteractionLagrangian(fermion_vertex, :u)

∂xψ₂ = partial(ψ₂, :x)
derivative_vertex = ψ₁ * ∂xψ₂ * bar(ψ₁) * bar(∂xψ₂)
L_p = InteractionLagrangian(derivative_vertex, :γp)

println(
    "workload\tstatistics\torder\tedges\tterms\tbaseline_states\tdirect_states\t",
    "baseline_completed\tdirect_completed\tbaseline_transitions\tdirect_transitions\t",
    "baseline_canon\tdirect_canon\tbaseline_pruned\tdirect_pruned\tonepi_pruned\toutputs\t",
    "baseline_s\tdirect_s\ttime_ratio\tbaseline_bytes\tdirect_bytes\talloc_ratio",
)

benchmark_onepi_workload("boson_g2", Boson, L_b, Val(2), Val(5); samples=3)
benchmark_onepi_workload("boson_g3", Boson, L_b, Val(3), Val(7); samples=2)
benchmark_onepi_workload("boson_g4", Boson, L_b, Val(4), Val(9); samples=1)
benchmark_onepi_workload("boson_gamma2", Boson, L_γ, Val(2), Val(5); samples=3)
benchmark_onepi_workload("boson_gamma3", Boson, L_γ, Val(3), Val(7); samples=2)
benchmark_onepi_workload(
    "fermion_u2", Fermion, L_f, Val(2), Val(5); simplify=false, samples=3
)
benchmark_onepi_workload(
    "fermion_u3", Fermion, L_f, Val(3), Val(7); simplify=false, samples=2
)
benchmark_onepi_workload(
    "fermion_derivative2", Fermion, L_p, Val(2), Val(5); simplify=false, samples=3
)
benchmark_onepi_workload(
    "fermion_derivative3", Fermion, L_p, Val(3), Val(7); simplify=false, samples=2
)
