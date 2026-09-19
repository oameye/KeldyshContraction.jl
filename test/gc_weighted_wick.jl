using KeldyshContraction, Test
import KeldyshContraction as KC

function signed_pairing_weights(pairings)
    weights = Dict{Any,Int}()
    for (pairing, topology, multiplicity) in pairings
        key = (pairing.contractions, topology)
        weight = Int(pairing.sign) * multiplicity
        weights[key] = get(weights, key, 0) + weight
    end
    filter!(pair -> !iszero(last(pair)), weights)
    return weights
end

@testset "GraphCombinations weighted Wick backend" begin
    @testset "bosonic exactness" begin
        @qfields gcwick_ϕ::Boson
        c, q = gcwick_ϕ[Classical], gcwick_ϕ[Quantum]
        elastic =
            -(0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2))

        expression = c(KC.Out()) * bar(q)(KC.In()) * elastic
        for term in KC.terms(expression)
            direct = KC._wick_contraction(
                term.args_nc, Val(3), Val(0); regularise=true, simplify=false
            )
            gc, stats = KC._gc_wick_contraction_with_stats(
                term.args_nc, Val(3), Val(0); regularise=true, simplify=false
            )
            @test signed_pairing_weights(gc) == signed_pairing_weights(direct)
            @test stats.canonicalization_calls >= 1
        end

        L = InteractionLagrangian(elastic)
        second_order = c(KC.Out()) * bar(c)(KC.In()) * L(1).lagrangian * L(2).lagrangian
        for term in KC.terms(second_order)
            direct = KC._wick_contraction(
                term.args_nc, Val(5), Val(1); regularise=true, simplify=true
            )
            gc, _ = KC._gc_wick_contraction_with_stats(
                term.args_nc, Val(5), Val(1); regularise=true, simplify=true
            )
            @test signed_pairing_weights(gc) == signed_pairing_weights(direct)
        end
    end

    @testset "fermionic cofactor cancellation" begin
        @qfields gcwick_pair_ψ::Fermion
        ψ₁ = gcwick_pair_ψ[One](KC.Bulk(1))
        ψ₂ = gcwick_pair_ψ[Two](KC.Bulk(1))
        bψ₂a = bar(gcwick_pair_ψ[Two])(KC.Bulk(2))
        bψ₂b = bar(gcwick_pair_ψ[Two])(KC.Bulk(3))
        term = ψ₁ * ψ₂ * bψ₂a * bψ₂b
        args = copy(KC.fields(term))

        direct = KC._wick_contraction(
            args, Val(2), Val(0); regularise=false, simplify=false
        )
        gc, stats = KC._gc_wick_contraction_with_stats(
            args, Val(2), Val(0); regularise=false, simplify=false
        )
        @test signed_pairing_weights(gc) == signed_pairing_weights(direct)
        @test isempty(gc)
        @test stats.automorphisms > 1
    end

    @testset "fermionic second-order relabeling transport" begin
        @qfields gcwick_ψ::Fermion
        ψ₁, ψ₂ = gcwick_ψ[One], gcwick_ψ[Two]
        vertex = ψ₁ * ψ₂ * bar(ψ₁) * bar(ψ₂)
        L = InteractionLagrangian(vertex)
        external_fields = KC.propagator_fields(L, nothing)
        products = KC.propagator_external_products(Fermion, external_fields...)

        saw_nontrivial_automorphism = false
        saw_merged_transition = false
        for in_out in products
            expression = in_out * L(1).lagrangian * L(2).lagrangian
            for term in KC.terms(expression)
                direct = KC._wick_contraction(
                    term.args_nc, Val(5), Val(1); regularise=false, simplify=false
                )
                gc, stats = KC._gc_wick_contraction_with_stats(
                    term.args_nc, Val(5), Val(1); regularise=false, simplify=false
                )
                @test signed_pairing_weights(gc) == signed_pairing_weights(direct)
                saw_nontrivial_automorphism |= stats.automorphisms > 1
                saw_merged_transition |= stats.merged_transitions > 0
            end
        end
        @test saw_nontrivial_automorphism
        @test saw_merged_transition
    end
end
