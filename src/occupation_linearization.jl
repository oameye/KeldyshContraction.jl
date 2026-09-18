"""
Canonical Fréchet derivative of an occupation polynomial.

Each entry maps one varied occupation atom `δn_a` to the exact background polynomial
multiplying that variation. Repeated factors are differentiated with their exact integer
multiplicity, and equal variation atoms are merged canonically.
"""
struct OccupationLinearization{C<:Number,S<:Statistics}
    terms::Vector{Pair{OccupationAtom{S},OccupationPolynomial{C,S}}}

    function OccupationLinearization{C,S}(
        terms::Vector{Pair{OccupationAtom{S},OccupationPolynomial{C,S}}}, ::Val{:raw}
    ) where {C<:Number,S<:Statistics}
        return new{C,S}(terms)
    end
end

function _canonical_occupation_linearization(
    terms::Vector{Pair{OccupationAtom{S},OccupationPolynomial{C,S}}}
) where {C<:Number,S<:Statistics}
    isempty(terms) && return OccupationLinearization{C,S}(terms, Val(:raw))
    sort!(terms; by=first)
    out = Pair{OccupationAtom{S},OccupationPolynomial{C,S}}[]
    sizehint!(out, length(terms))
    for (atom, polynomial) in terms
        iszero(polynomial) && continue
        if !isempty(out) && isequal(first(out[end]), atom)
            combined = last(out[end]) + polynomial
            pop!(out)
            iszero(combined) || push!(out, atom => combined)
        else
            push!(out, atom => polynomial)
        end
    end
    return OccupationLinearization{C,S}(out, Val(:raw))
end

function OccupationLinearization{C,S}(
    terms::Vector{Pair{OccupationAtom{S},OccupationPolynomial{C,S}}}
) where {C<:Number,S<:Statistics}
    return _canonical_occupation_linearization(copy(terms))
end

Base.length(linearization::OccupationLinearization) = length(linearization.terms)
Base.isempty(linearization::OccupationLinearization) = isempty(linearization.terms)
Base.iszero(linearization::OccupationLinearization) = isempty(linearization.terms)
Base.iterate(linearization::OccupationLinearization) = iterate(linearization.terms)
function Base.iterate(linearization::OccupationLinearization, state)
    return iterate(linearization.terms, state)
end
function Base.eltype(::Type{OccupationLinearization{C,S}}) where {C,S}
    return Pair{OccupationAtom{S},OccupationPolynomial{C,S}}
end
Base.IteratorSize(::Type{<:OccupationLinearization}) = Base.HasLength()

"""Return the canonical variation-atom to background-polynomial entries."""
occupation_linearization_terms(linearization::OccupationLinearization) = linearization.terms

"""
    occupation_linearization(polynomial)

Construct the exact first Fréchet derivative of a canonical occupation polynomial.
"""
function occupation_linearization(
    polynomial::OccupationPolynomial{C,S}
) where {C<:Number,S<:Statistics}
    D = promote_type(C, Int)
    raw = Pair{OccupationAtom{S},OccupationPolynomial{D,S}}[]

    for (monomial, coefficient) in polynomial
        factors = monomial.factors
        i = firstindex(factors)
        while i <= lastindex(factors)
            atom = factors[i]
            j = i
            while j < lastindex(factors) && isequal(factors[j + 1], atom)
                j += 1
            end
            multiplicity = j - i + 1
            residual_factors = copy(factors)
            deleteat!(residual_factors, i)
            residual = OccupationMonomial(residual_factors)
            derivative_coefficient = convert(D, coefficient) * convert(D, multiplicity)
            contribution = OccupationPolynomial{D,S}([residual => derivative_coefficient])
            push!(raw, atom => contribution)
            i = j + 1
        end
    end

    return _canonical_occupation_linearization(raw)
end

"""
Occupation-space response coefficients evaluated at a supplied background state.

Each entry retains the exact variation channel `δn_a` while replacing its background occupation
polynomial by the corresponding value at `n̄`.
"""
struct BackgroundOccupationLinearization{C<:Number,S<:Statistics}
    terms::Vector{Pair{OccupationAtom{S},C}}

    function BackgroundOccupationLinearization{C,S}(
        terms::Vector{Pair{OccupationAtom{S},C}}, ::Val{:raw}
    ) where {C<:Number,S<:Statistics}
        return new{C,S}(terms)
    end
end

function _canonical_background_occupation_linearization(
    terms::Vector{Pair{OccupationAtom{S},C}}
) where {C<:Number,S<:Statistics}
    isempty(terms) && return BackgroundOccupationLinearization{C,S}(terms, Val(:raw))
    sort!(terms; by=first)
    out = Pair{OccupationAtom{S},C}[]
    sizehint!(out, length(terms))
    for (atom, coefficient) in terms
        iszero(coefficient) && continue
        if !isempty(out) && isequal(first(out[end]), atom)
            combined = last(out[end]) + coefficient
            pop!(out)
            iszero(combined) || push!(out, atom => combined)
        else
            push!(out, atom => coefficient)
        end
    end
    return BackgroundOccupationLinearization{C,S}(out, Val(:raw))
end

function BackgroundOccupationLinearization{C,S}(
    terms::Vector{Pair{OccupationAtom{S},C}}
) where {C<:Number,S<:Statistics}
    return _canonical_background_occupation_linearization(copy(terms))
end

Base.length(linearization::BackgroundOccupationLinearization) = length(linearization.terms)
function Base.isempty(linearization::BackgroundOccupationLinearization)
    return isempty(linearization.terms)
end
Base.iszero(linearization::BackgroundOccupationLinearization) = isempty(linearization.terms)
function Base.iterate(linearization::BackgroundOccupationLinearization)
    return iterate(linearization.terms)
end
function Base.iterate(linearization::BackgroundOccupationLinearization, state)
    return iterate(linearization.terms, state)
end
function Base.eltype(::Type{BackgroundOccupationLinearization{C,S}}) where {C,S}
    return Pair{OccupationAtom{S},C}
end
Base.IteratorSize(::Type{<:BackgroundOccupationLinearization}) = Base.HasLength()

"""Return the exact variation-channel coefficients evaluated at the background state."""
function background_occupation_linearization_terms(
    linearization::BackgroundOccupationLinearization
)
    return linearization.terms
end

"""Evaluate an exact occupation-space Fréchet derivative at a supplied background state."""
function evaluate_occupation_linearization(
    linearization::OccupationLinearization{C,S}, background::OccupationBackground{V}
) where {C<:Number,S<:Statistics,V<:Number}
    D = promote_type(C, V)
    terms = Pair{OccupationAtom{S},D}[]
    sizehint!(terms, length(linearization))
    for (atom, polynomial) in linearization
        coefficient = convert(D, evaluate_occupation_polynomial(polynomial, background))
        push!(terms, atom => coefficient)
    end
    return BackgroundOccupationLinearization{D,S}(terms)
end
