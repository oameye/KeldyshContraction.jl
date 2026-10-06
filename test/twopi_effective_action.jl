using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields twopi_boson::Boson

@testset "formal 2PI effective action" begin
    c = twopi_boson[Classical]
    q = twopi_boson[Quantum]
    vertex = c * q * bar(c) * bar(q)
    L = InteractionLagrangian(vertex, :λ)

    Γ2 = @inferred TwoPIEffectiveAction(L, Val(1), Val(2))
    @test KC.order(Γ2) == 1
    @test KC.statistics(Γ2) === Boson
    @test KC.parameters(Γ2) == KC.parameter_monomial(:λ)
    @test KC.field_families(Γ2) == KC.field_families(L)
    @test !iszero(Γ2)

    terms = KC.twopi_terms(Γ2)
    @test length(terms) == 2
    @test all(value == convert(KC.ComplexRationals, -1) for value in values(terms))
    @test all(
        KC.is_two_particle_irreducible(collect(KC.twopi_contractions(diagram))) for
        diagram in keys(terms)
    )
    @test all(
        all(
            KC.is_bulk,
            (
                field for contraction in KC.twopi_contractions(diagram) for
                field in contraction
            ),
        ) for diagram in keys(terms)
    )

    # The physical perturbative vacuum path removes both formal contributions: one contains Gqq,
    # while the other is a causal RA zero loop. Γ₂ must retain them before functional
    # differentiation because a vanishing physical component may have a nonzero derivative.
    physical_vacuum = KC._wick_contraction(vertex, Val(2); simplify=false)
    KC.filter_nonzero!(physical_vacuum)
    @test iszero(physical_vacuum)

    formal_contractions = [
        contraction for diagram in keys(terms) for
        contraction in KC.twopi_contractions(diagram)
    ]
    @test any(KC.is_qq_contraction, formal_contractions)
end

@qfields twopi_fermion::Fermion

@testset "fermionic 2PI Wick parity" begin
    c = twopi_fermion[Classical]
    q = twopi_fermion[Quantum]
    vertex = c * q * bar(c) * bar(q)
    L = InteractionLagrangian(vertex, :u)

    Γ2 = @inferred TwoPIEffectiveAction(L, Val(1), Val(2))
    @test KC.statistics(Γ2) === Fermion
    @test KC.parameters(Γ2) == KC.parameter_monomial(:u)

    coefficients = collect(values(KC.twopi_terms(Γ2)))
    @test length(coefficients) == 2
    @test coefficients[1] == -coefficients[2]
    @test all(value -> !iszero(value), coefficients)
end

@testset "2PI static shape validation" begin
    c = twopi_boson[Classical]
    q = twopi_boson[Quantum]
    L = InteractionLagrangian(c * q * bar(c) * bar(q), :λ)

    @test_throws ArgumentError TwoPIEffectiveAction(L, Val(0), Val(0))
    @test_throws ArgumentError TwoPIEffectiveAction(L, Val(1), Val(3))
end
