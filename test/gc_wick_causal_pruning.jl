using KeldyshContraction, Test
import KeldyshContraction as KC

function causal_pairing_weights(pairings)
    weights = Dict{Any,Int}()
    for (pairing, topology, multiplicity) in pairings
        key = (pairing.contractions, topology)
        weight = Int(pairing.sign) * multiplicity
        weights[key] = get(weights, key, 0) + weight
    end
    filter!(pair -> !iszero(last(pair)), weights)
    return weights
end

function certify_causal_pruning(args_nc, ::Val{E}, ::Val{E2}; kwargs...) where {E,E2}
    direct = KC._wick_contraction(args_nc, Val(E), Val(E2); kwargs...)
    unpruned, unpruned_stats = KC._gc_wick_contraction_with_stats(
        args_nc, Val(E), Val(E2); causal_pruning=false, kwargs...
    )
    pruned, pruned_stats = KC._gc_wick_contraction_with_stats(
        args_nc, Val(E), Val(E2); causal_pruning=true, kwargs...
    )
    expected = causal_pairing_weights(direct)
    @test causal_pairing_weights(unpruned) == expected
    @test causal_pairing_weights(pruned) == expected
    return unpruned_stats, pruned_stats
end

@testset "GraphCombinations causal Wick pruning" begin
    @qfields causal_ϕ::Boson
    c, q = causal_ϕ[Classical], causal_ϕ[Quantum]
    boson_vertex = -(
        0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2)
    )
    L_b = InteractionLagrangian(boson_vertex)
    external_b = c(KC.Out()) * bar(c)(KC.In())

    total_unpruned_canon = 0
    total_pruned_canon = 0
    total_pruned_transitions = 0
    for order in (3, 4)
        expression = external_b * prod(L_b(i).lagrangian for i in 1:order)
        E = 2order + 1
        E2 = KC.max_edges(order)
        for term in KC.terms(expression)
            unpruned_stats, pruned_stats = certify_causal_pruning(
                term.args_nc,
                Val(E),
                Val(E2);
                regularise=true,
                _set_reg_to_zero=false,
                simplify=true,
            )
            total_unpruned_canon += unpruned_stats.canonicalization_calls
            total_pruned_canon += pruned_stats.canonicalization_calls
            total_pruned_transitions += pruned_stats.pruned_transitions
        end
    end
    @test total_pruned_transitions > 0
    @test total_pruned_canon < total_unpruned_canon

    @qfields causal_ψ::Fermion
    ψ₁, ψ₂ = causal_ψ[One], causal_ψ[Two]
    fermion_vertex = ψ₁ * ψ₂ * bar(ψ₁) * bar(ψ₂)
    L_f = InteractionLagrangian(fermion_vertex)
    products = KC.propagator_external_products(
        Fermion, KC.propagator_fields(L_f, nothing)...
    )

    fermion_pruned = 0
    for in_out in products
        expression = in_out * L_f(1).lagrangian * L_f(2).lagrangian * L_f(3).lagrangian
        for term in KC.terms(expression)
            _, pruned_stats = certify_causal_pruning(
                term.args_nc,
                Val(7),
                Val(3);
                regularise=false,
                _set_reg_to_zero=false,
                simplify=false,
            )
            fermion_pruned += pruned_stats.pruned_transitions
        end
    end
    @test fermion_pruned > 0
end
