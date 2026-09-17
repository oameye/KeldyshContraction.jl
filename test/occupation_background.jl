using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields occupation_background_ϕ::Boson

@testset "occupation background public boundary" begin
    for name in (
        :OccupationBackground, :occupation_background_value, :evaluate_occupation_polynomial
    )
        VERSION >= v"1.11" && @test Base.ispublic(KC, name)
        VERSION >= v"1.11" && @test !Base.isexported(KC, name)
    end
end

@testset "typed occupation background evaluates exact polynomials" begin
    C = KC.ComplexRationals
    basis = KC.MomentumBasis(2)
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    nk = KC.OccupationPolynomial(KC.OccupationAtom(occupation_background_ϕ, k), one(C))
    nq = KC.OccupationPolynomial(KC.OccupationAtom(occupation_background_ϕ, q), one(C))
    k_atom = KC.OccupationAtom(occupation_background_ϕ, k)
    q_atom = KC.OccupationAtom(occupation_background_ϕ, q)

    background = KC.OccupationBackground(Rational{Int64}) do atom
        return isequal(atom, k_atom) ? 2 // 1 : 5 // 1
    end

    @test @inferred(KC.occupation_background_value(background, k_atom)) == 2 // 1
    @test @inferred(KC.occupation_background_value(background, q_atom)) == 5 // 1

    explicit_background = KC.OccupationBackground(_ -> 3 // 1, Rational{Int64})
    @test @inferred(KC.occupation_background_value(explicit_background, k_atom)) == 3 // 1

    polynomial = 2 * nk * nk * nq + 3 * nq
    @test @inferred(KC.evaluate_occupation_polynomial(polynomial, background)) == 55 // 1
end
