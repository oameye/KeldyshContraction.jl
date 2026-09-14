"""
Exact affine causal denominator used by canonical collision-frequency reduction.

The denominator represents

    a⋅ω + E + iη0

with exact rational loop-frequency coefficients, canonical symbolic `EnergyForm` `E`, and an
exact causal prescription `η`. Equivalent nonzero rational rescalings are canonicalized at the
term boundary, where the inverse scale is transferred into the term coefficient.
"""
struct CausalFrequencyDenominator{S<:Statistics}
    loop_coefficients::Vector{EnergyCoefficient}
    energy::EnergyForm{S}
    infinitesimal::EnergyCoefficient
end

function CausalFrequencyDenominator(
    loop_coefficients::AbstractVector{<:Rational},
    energy::EnergyForm{S},
    infinitesimal::Rational,
) where {S<:Statistics}
    return CausalFrequencyDenominator{S}(
        EnergyCoefficient[convert(EnergyCoefficient, c) for c in loop_coefficients],
        energy,
        convert(EnergyCoefficient, infinitesimal),
    )
end

function Base.isequal(
    a::CausalFrequencyDenominator{S}, b::CausalFrequencyDenominator{S}
) where {S<:Statistics}
    return isequal(a.loop_coefficients, b.loop_coefficients) &&
           isequal(a.energy, b.energy) &&
           a.infinitesimal == b.infinitesimal
end
Base.:(==)(a::CausalFrequencyDenominator, b::CausalFrequencyDenominator) = isequal(a, b)
function Base.hash(denominator::CausalFrequencyDenominator, h::UInt)
    return hash(
        denominator.infinitesimal,
        hash(
            denominator.energy,
            hash(denominator.loop_coefficients, hash(CausalFrequencyDenominator, h)),
        ),
    )
end

function Base.isless(
    a::CausalFrequencyDenominator{S}, b::CausalFrequencyDenominator{S}
) where {S<:Statistics}
    length(a.loop_coefficients) == length(b.loop_coefficients) ||
        return length(a.loop_coefficients) < length(b.loop_coefficients)
    @inbounds for i in eachindex(a.loop_coefficients, b.loop_coefficients)
        a.loop_coefficients[i] == b.loop_coefficients[i] ||
            return a.loop_coefficients[i] < b.loop_coefficients[i]
    end
    isequal(a.energy, b.energy) || return isless(a.energy, b.energy)
    return a.infinitesimal < b.infinitesimal
end

function _causal_frequency_scale(denominator::CausalFrequencyDenominator)
    for coefficient in denominator.loop_coefficients
        iszero(coefficient) || return coefficient
    end
    for (_, coefficient) in denominator.energy
        iszero(coefficient) || return coefficient
    end
    iszero(denominator.infinitesimal) || return denominator.infinitesimal
    return throw(
        ArgumentError("identically zero causal denominator has no canonical scale")
    )
end

"""
Return `(canonical, factor)` for one affine causal denominator.

If the supplied denominator is `c D₀`, the canonical denominator has first nonzero exact
coefficient one and `factor == 1/c`, so `1/(c D₀) = factor / D₀`. The infinitesimal is scaled
with the same rational factor as the real affine form; its orientation therefore remains part of
the canonical denominator identity.
"""
function canonical_causal_frequency_denominator(
    denominator::CausalFrequencyDenominator{S}
) where {S<:Statistics}
    scale = _causal_frequency_scale(denominator)
    inverse_scale = inv(scale)
    canonical = CausalFrequencyDenominator{S}(
        EnergyCoefficient[c * inverse_scale for c in denominator.loop_coefficients],
        inverse_scale * denominator.energy,
        denominator.infinitesimal * inverse_scale,
    )
    return canonical, inverse_scale
end

"""One canonical product of affine causal frequency denominators."""
struct CausalFrequencyTerm{C<:Number,S<:Statistics}
    coefficient::C
    denominators::Vector{CausalFrequencyDenominator{S}}

    function CausalFrequencyTerm{C,S}(
        coefficient::C,
        denominators::Vector{CausalFrequencyDenominator{S}},
        ::Val{:canonical},
    ) where {C<:Number,S<:Statistics}
        return new{C,S}(coefficient, denominators)
    end
end

function CausalFrequencyTerm(
    coefficient::C, denominators::AbstractVector{CausalFrequencyDenominator{S}}
) where {C<:Number,S<:Statistics}
    D = promote_type(C, EnergyCoefficient)
    canonical_coefficient = convert(D, coefficient)
    canonical = CausalFrequencyDenominator{S}[]
    sizehint!(canonical, length(denominators))
    nfrequencies =
        isempty(denominators) ? nothing : length(first(denominators).loop_coefficients)
    basis_size =
        isempty(denominators) ? nothing : energy_basis_size(first(denominators).energy)

    for denominator in denominators
        if nfrequencies !== nothing
            length(denominator.loop_coefficients) == nfrequencies || throw(
                DimensionMismatch("causal denominators use different loop-frequency bases"),
            )
            energy_basis_size(denominator.energy) == basis_size ||
                throw(DimensionMismatch("causal denominators use different energy bases"))
        end
        normalized, factor = canonical_causal_frequency_denominator(denominator)
        canonical_coefficient *= convert(D, factor)
        push!(canonical, normalized)
    end
    sort!(canonical)
    return CausalFrequencyTerm{D,S}(canonical_coefficient, canonical, Val(:canonical))
end

causal_frequency_coefficient(term::CausalFrequencyTerm) = term.coefficient
causal_frequency_denominators(term::CausalFrequencyTerm) = term.denominators

function Base.isequal(a::CausalFrequencyTerm{C,S}, b::CausalFrequencyTerm{C,S}) where {C,S}
    return a.coefficient == b.coefficient && isequal(a.denominators, b.denominators)
end
Base.:(==)(a::CausalFrequencyTerm, b::CausalFrequencyTerm) = isequal(a, b)
function Base.hash(term::CausalFrequencyTerm, h::UInt)
    return hash(term.denominators, hash(term.coefficient, hash(CausalFrequencyTerm, h)))
end

function _causal_denominator_vectors_isless(
    a::Vector{CausalFrequencyDenominator{S}}, b::Vector{CausalFrequencyDenominator{S}}
) where {S<:Statistics}
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        isequal(a[i], b[i]) || return isless(a[i], b[i])
    end
    return length(a) < length(b)
end

"""
Canonical sum of causal frequency terms.

Terms with identical canonical denominator multisets are merged before any contour decision is
made. This is the semantic boundary required for cancellations that exist only in the complete
collision expression rather than in an isolated spectral/dispersive term.
"""
struct CausalFrequencyExpression{C<:Number,S<:Statistics}
    terms::Vector{CausalFrequencyTerm{C,S}}

    function CausalFrequencyExpression{C,S}(
        terms::Vector{CausalFrequencyTerm{C,S}}, ::Val{:canonical}
    ) where {C<:Number,S<:Statistics}
        return new{C,S}(terms)
    end
end

function CausalFrequencyExpression(
    terms::AbstractVector{CausalFrequencyTerm{C,S}}
) where {C<:Number,S<:Statistics}
    ordered = collect(terms)
    sort!(
        ordered;
        lt=(a, b) -> _causal_denominator_vectors_isless(a.denominators, b.denominators),
    )

    merged = CausalFrequencyTerm{C,S}[]
    sizehint!(merged, length(ordered))
    for term in ordered
        iszero(term.coefficient) && continue
        if !isempty(merged) && isequal(merged[end].denominators, term.denominators)
            coefficient = merged[end].coefficient + term.coefficient
            denominators = merged[end].denominators
            pop!(merged)
            iszero(coefficient) || push!(
                merged,
                CausalFrequencyTerm{C,S}(coefficient, denominators, Val(:canonical)),
            )
        else
            push!(merged, term)
        end
    end
    return CausalFrequencyExpression{C,S}(merged, Val(:canonical))
end

causal_frequency_terms(expression::CausalFrequencyExpression) = expression.terms
Base.isempty(expression::CausalFrequencyExpression) = isempty(expression.terms)
Base.iszero(expression::CausalFrequencyExpression) = isempty(expression)

function Base.isequal(
    a::CausalFrequencyExpression{C,S}, b::CausalFrequencyExpression{C,S}
) where {C,S}
    return isequal(a.terms, b.terms)
end
Base.:(==)(a::CausalFrequencyExpression, b::CausalFrequencyExpression) = isequal(a, b)
function Base.hash(expression::CausalFrequencyExpression, h::UInt)
    return hash(expression.terms, hash(CausalFrequencyExpression, h))
end

@inline function _frequency_incidence(term::CausalFrequencyTerm, frequency_index::Int)::Int
    return count(
        denominator -> !iszero(denominator.loop_coefficients[frequency_index]),
        term.denominators,
    )
end

function _check_causal_frequency_index(
    expression::CausalFrequencyExpression, frequency_index::Int
)::Nothing
    isempty(expression) && return nothing
    denominators = first(expression.terms).denominators
    isempty(denominators) &&
        throw(ArgumentError("frequency expression has no denominators"))
    nfrequencies = length(first(denominators).loop_coefficients)
    checkbounds(1:nfrequencies, frequency_index)
    return nothing
end

function _leading_frequency_term(
    term::CausalFrequencyTerm{C,S}, frequency_index::Int
) where {C<:Number,S<:Statistics}
    D = promote_type(C, EnergyCoefficient)
    coefficient = convert(D, term.coefficient)
    spectators = CausalFrequencyDenominator{S}[]
    sizehint!(spectators, length(term.denominators))
    for denominator in term.denominators
        frequency_coefficient = denominator.loop_coefficients[frequency_index]
        if iszero(frequency_coefficient)
            push!(spectators, denominator)
        else
            coefficient /= convert(D, frequency_coefficient)
        end
    end
    return CausalFrequencyTerm(coefficient, spectators)
end

"""
Return the first exact large-`ω_j` asymptotic coefficient of a complete causal expression.

The integer is the smallest denominator incidence `p` among terms, corresponding to a candidate
`ω_j^-p` tail. The returned expression is the exact coefficient of that power as a rational
function of all spectator frequencies. It is canonicalized and combined across the complete
input expression, so a vanishing result proves cancellation of that tail globally.
"""
function leading_frequency_asymptotic(
    expression::CausalFrequencyExpression{C,S}, frequency_index::Int
) where {C<:Number,S<:Statistics}
    _check_causal_frequency_index(expression, frequency_index)
    D = promote_type(C, EnergyCoefficient)
    if isempty(expression)
        empty_terms = CausalFrequencyTerm{D,S}[]
        return typemax(Int), CausalFrequencyExpression{D,S}(empty_terms, Val(:canonical))
    end

    order = typemax(Int)
    for term in expression.terms
        order = min(order, _frequency_incidence(term, frequency_index))
    end

    leading_terms = CausalFrequencyTerm{D,S}[]
    for term in expression.terms
        _frequency_incidence(term, frequency_index) == order || continue
        push!(leading_terms, _leading_frequency_term(term, frequency_index))
    end
    return order, CausalFrequencyExpression(leading_terms)
end

"""
Return a certified lower bound on decay in one loop frequency.

If the leading `ω^-p` coefficient cancels in the complete canonical expression, the result is at
least `p+1`. Only the first cancellation is needed to distinguish the marginal `1/ω` case from
the ordinary contour-safe `O(1/ω²)` case; deeper asymptotic expansion can be added without
changing the IR.
"""
function frequency_decay_lower_bound(
    expression::CausalFrequencyExpression{C,S}, frequency_index::Int
)::Int where {C<:Number,S<:Statistics}
    isempty(expression) && return typemax(Int)
    order, leading = leading_frequency_asymptotic(expression, frequency_index)
    return isempty(leading) ? order + 1 : order
end

"""Return whether the complete expression is certified to have a vanishing large semicircle."""
function causal_contour_safe_at_infinity(
    expression::CausalFrequencyExpression{C,S}, frequency_index::Int
)::Bool where {C<:Number,S<:Statistics}
    return frequency_decay_lower_bound(expression, frequency_index) >= 2
end
