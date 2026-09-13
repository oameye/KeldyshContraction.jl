using KeldyshContraction, Test
import KeldyshContraction as KC

function signed_pairing_weights(pairings)
    weights = Dict()
    for (pairing, topology, multiplicity) in pairings
        key = (pairing.contractions, topology)
        signed = Int(pairing.sign) * multiplicity
        weights[key] = get(weights, key, 0) + signed
    end
    filter!(pair -> !iszero(last(pair)), weights)
    return weights
end

@testset "GraphCombinations bosonic Wick backend" begin
    @qfields ϕ::Boson
    c, q = ϕ[Classical], ϕ[Quantum]
    elastic = -(0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2))

    @testset "first-order physical term" begin
        expression = c(Out()) * bar(q)(In()) * elastic
        for term in KC.terms(expression)
            direct = KC._wick_contraction(
                term.args_nc, Val(3), Val(0); regularise=true, simplify=false
            )
            gc, stats = KC._gc_bosonic_wick_contraction_with_stats(
                term.args_nc, Val(3), Val(0); regularise=true, simplify=false
            )
            @test signed_pairing_weights(gc) == signed_pairing_weights(direct)
            @test stats.canonicalization_calls >= 1
        end
    end

    @testset "second-order interaction products" begin
        L = InteractionLagrangian(elastic)
        expression = c(Out()) * bar(c)(In()) * L(1).lagrangian * L(2).lagrangian
        for term in KC.terms(expression)
            direct = KC._wick_contraction(
                term.args_nc, Val(5), Val(1); regularise=true, simplify=true
            )
            gc = KC._gc_bosonic_wick_contraction(
                term.args_nc, Val(5), Val(1); regularise=true, simplify=true
            )
            @test signed_pairing_weights(gc) == signed_pairing_weights(direct)
        end
    end
end
