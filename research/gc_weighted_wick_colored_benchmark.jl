using KeldyshContraction
using Combinatorics
import KeldyshContraction as KC

function colored_semantic_weights(pairings)
    weights = Dict{Any,Int}()
    for (pairing, topology, multiplicity) in pairings
        key = (pairing.contractions, topology)
        weight = Int(pairing.sign) * multiplicity
        weights[key] = get(weights, key, 0) + weight
    end
    filter!(pair -> !iszero(last(pair)), weights)
    return weights
end

function colored_workload_terms(in_out, L, ::Val{O}) where {O}
    result = Vector{typeof(copy(in_out.args_nc))}()
    for coefficients in Combinatorics.multiexponents(length(L.lagrangian), O)
        idxs = KC.indices_from_counts(coefficients)
        mult = Combinatorics.multinomial(coefficients...)
        qmul = mult * prod(L(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))
        push!(result, copy((in_out * qmul).args_nc))
    end
    return result
end

function colored_interaction_terms(::Type{S}, L, ::Val{O}) where {S<:KC.Statistics,O}
    products = KC.propagator_external_products(S, KC.propagator_fields(L, nothing)...)
    result = Vector{typeof(copy(first(products).args_nc))}()
    for in_out in products
        append!(result, colored_workload_terms(in_out, L, Val(O)))
    end
    return result
end

function colored_direct_run(
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

function colored_gc_run(
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

function colored_best_measurement(f; samples=5)
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

function benchmark_colored_workload(
    name, ::Type{S}, L, ::Val{O}, ::Val{E}; simplify=true, samples=5
) where {S,O,E}
    terms = colored_interaction_terms(S, L, Val(O))
    E2 = KC.max_edges(O)
    regularise = KC.should_regularise(L.lagrangian)
    kwargs = (; regularise, _set_reg_to_zero=true, simplify)

    for args_nc in terms
        direct = KC._wick_contraction(args_nc, Val(E), Val(E2); kwargs...)
        gc = KC._gc_wick_contraction(args_nc, Val(E), Val(E2); kwargs...)
        colored_semantic_weights(gc) == colored_semantic_weights(direct) ||
            error("GC/direct semantic mismatch in colored crossover benchmark: $name")
    end

    direct = () -> colored_direct_run(terms, Val(E), Val(E2); kwargs...)
    gc = () -> colored_gc_run(terms, Val(E), Val(E2); kwargs...)
    direct_outputs = direct()
    gc_outputs = gc()
    direct_outputs == gc_outputs || error("GC/direct output-count mismatch: $name")

    direct_time, direct_bytes = colored_best_measurement(direct; samples)
    gc_time, gc_bytes = colored_best_measurement(gc; samples)
    println(
        join(
            (
                name,
                string(S),
                O,
                E,
                length(terms),
                direct_outputs,
                direct_time,
                gc_time,
                gc_time / direct_time,
                direct_bytes,
                gc_bytes,
                gc_bytes / direct_bytes,
            ),
            '\t',
        ),
    )
    return nothing
end

@qfields colored_ϕ::Boson
c, q = colored_ϕ[Classical], colored_ϕ[Quantum]
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

@qfields colored_ψ::Fermion
ψ₁, ψ₂ = colored_ψ[One], colored_ψ[Two]
∂xψ₂ = partial(ψ₂, :x)
derivative_vertex = ψ₁ * ∂xψ₂ * bar(ψ₁) * bar(∂xψ₂)
L_p = InteractionLagrangian(derivative_vertex, :γp)

println(
    "workload\tstatistics\torder\tedges\tterms\toutputs\tdirect_s\tgc_s\ttime_ratio\t",
    "direct_bytes\tgc_bytes\talloc_ratio",
)
benchmark_colored_workload("boson_gamma2", Boson, L_γ, Val(2), Val(5); samples=5)
benchmark_colored_workload("boson_gamma3", Boson, L_γ, Val(3), Val(7); samples=3)
benchmark_colored_workload(
    "fermion_derivative2", Fermion, L_p, Val(2), Val(5); simplify=false, samples=5
)
benchmark_colored_workload(
    "fermion_derivative3", Fermion, L_p, Val(3), Val(7); simplify=false, samples=5
)
