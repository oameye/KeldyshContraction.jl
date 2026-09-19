using KeldyshContraction, Test
import KeldyshContraction as KC

@testset "GraphCombinations Wick production crossover" begin
    @test KC._use_gc_wick(Boson, Val(5))
    @test KC._use_gc_wick(Boson, Val(7))
    @test !KC._use_gc_wick(Fermion, Val(5))
    @test KC._use_gc_wick(Fermion, Val(7))

    @qfields routing_ϕ::Boson
    c, q = routing_ϕ[Classical], routing_ϕ[Quantum]
    boson_vertex = -(
        0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2)
    )
    L_b = InteractionLagrangian(boson_vertex)
    boson_term = first(KC.terms(c(KC.Out()) * bar(c)(KC.In()) * L_b(1).lagrangian * L_b(2).lagrangian))
    direct_b = KC._wick_contraction(
        boson_term.args_nc, Val(5), Val(1); regularise=true, simplify=true
    )
    routed_b = KC._production_wick_pairings(
        boson_term.args_nc, Val(5), Val(1); regularise=true, simplify=true
    )
    @test signed_pairing_weights(routed_b) == signed_pairing_weights(direct_b)

    @qfields routing_ψ::Fermion
    ψ₁, ψ₂ = routing_ψ[One], routing_ψ[Two]
    fermion_vertex = ψ₁ * ψ₂ * bar(ψ₁) * bar(ψ₂)
    L_f = InteractionLagrangian(fermion_vertex)
    products = KC.propagator_external_products(Fermion, KC.propagator_fields(L_f, nothing)...)
    fermion_term2 = first(KC.terms(first(products) * L_f(1).lagrangian * L_f(2).lagrangian))
    direct_f2 = KC._wick_contraction(
        fermion_term2.args_nc, Val(5), Val(1); regularise=false, simplify=false
    )
    routed_f2 = KC._production_wick_pairings(
        fermion_term2.args_nc, Val(5), Val(1); regularise=false, simplify=false
    )
    @test signed_pairing_weights(routed_f2) == signed_pairing_weights(direct_f2)

    fermion_term3 = first(
        KC.terms(first(products) * L_f(1).lagrangian * L_f(2).lagrangian * L_f(3).lagrangian)
    )
    direct_f3 = KC._wick_contraction(
        fermion_term3.args_nc, Val(7), Val(3); regularise=false, simplify=false
    )
    routed_f3 = KC._production_wick_pairings(
        fermion_term3.args_nc, Val(7), Val(3); regularise=false, simplify=false
    )
    @test signed_pairing_weights(routed_f3) == signed_pairing_weights(direct_f3)
end
