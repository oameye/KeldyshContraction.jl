"""
Typed occupation background used to evaluate exact occupation polynomials.

The evaluator receives an `OccupationAtom` and must return a value convertible to `V`. Construct
as `OccupationBackground(f, V)`; the evaluator-first order also supports Julia `do` blocks while
keeping the scalar value type explicit for inference.
"""
struct OccupationBackground{V<:Number,F}
    evaluator::F
end

function OccupationBackground(evaluator::F, ::Type{V}) where {V<:Number,F}
    return OccupationBackground{V,F}(evaluator)
end

"""Evaluate one occupation atom in a typed supplied background state."""
function occupation_background_value(
    background::OccupationBackground{V}, atom::OccupationAtom
)::V where {V<:Number}
    return convert(V, background.evaluator(atom))
end

"""Evaluate an exact occupation polynomial in a supplied typed background state."""
function evaluate_occupation_polynomial(
    polynomial::OccupationPolynomial{C,S}, background::OccupationBackground{V}
) where {C<:Number,S<:Statistics,V<:Number}
    D = promote_type(C, V)
    value = zero(D)
    for (monomial, coefficient) in polynomial
        term = convert(D, coefficient)
        for atom in monomial
            term *= convert(D, occupation_background_value(background, atom))
        end
        value += term
    end
    return value
end
