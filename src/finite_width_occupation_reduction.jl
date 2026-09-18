"""
Resolved finite-width collision after exact statistical-to-occupation lowering.

`terms` contains only analytically evaluated finite-width frequency sectors. The original
`FiniteWidthFrequencyTerm` remains the key so the exact spectral-reduction provenance and routed
kinematics stay inspectable. Unsupported offset and external-distribution sectors remain separate
and are not assigned an occupation-space value.

The supplied spectral model has already been evaluated before this stage and is therefore treated
as fixed. In particular, this representation does not yet include the response of a
state-dependent linewidth `Γ[n]`.
"""
struct FiniteWidthOccupationCollision{
    C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext
}
    terms::Dict{FiniteWidthFrequencyTerm{S,E1,E2},OccupationPolynomial{C,S}}
    unsupported_offset::Dict{FiniteWidthFrequencyTerm{S,E1,E2},C}
    unsupported_distribution::Dict{FiniteWidthFrequencyTerm{S,E1,E2},C}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::FiniteWidthOccupationCollision{C,S,O}) where {C,S,O} = O
statistics(::FiniteWidthOccupationCollision{C,S}) where {C,S} = S
parameters(collision::FiniteWidthOccupationCollision) = collision.parameter
function gradient_order(
    ::FiniteWidthOccupationCollision{C,S,O,E1,E2,G}
) where {C,S,O,E1,E2,G}
    return Val(G)
end
wigner_context(collision::FiniteWidthOccupationCollision) = collision.context
target_family(collision::FiniteWidthOccupationCollision) = collision.target
finite_width_occupation_terms(collision::FiniteWidthOccupationCollision) = collision.terms
function finite_width_unsupported_offset_terms(collision::FiniteWidthOccupationCollision)
    return collision.unsupported_offset
end
function finite_width_unsupported_distribution_terms(
    collision::FiniteWidthOccupationCollision
)
    return collision.unsupported_distribution
end

function Base.isempty(collision::FiniteWidthOccupationCollision)
    return isempty(collision.terms) &&
           isempty(collision.unsupported_offset) &&
           isempty(collision.unsupported_distribution)
end

function _push_finite_width_occupation!(
    out::Dict{FiniteWidthFrequencyTerm{S,E1,E2},OccupationPolynomial{D,S}},
    term::FiniteWidthFrequencyTerm{S,E1,E2},
    occupation::OccupationPolynomial{D,S},
) where {D<:Number,S<:Statistics,E1,E2}
    iszero(occupation) && return out
    if haskey(out, term)
        combined = out[term] + occupation
        if iszero(combined)
            delete!(out, term)
        else
            out[term] = combined
        end
    else
        out[term] = occupation
    end
    return out
end

function _append_finite_width_occupation_component!(
    out::Dict{FiniteWidthFrequencyTerm{S,E1,E2},OccupationPolynomial{D,S}},
    expression::FiniteWidthEvaluatedExpression{C,S,E1,E2,G,Ctx},
    target::FieldFamily{S},
    include_external::Bool,
) where {C<:Number,D<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext}
    for (term, coefficient) in finite_width_resolved_terms(expression)
        statistical = _collision_statistical_monomial(
            finite_width_source_term(term), target, include_external
        )
        polynomial = StatisticalPolynomial{C,S}([statistical => coefficient])
        occupation = occupation_collision_polynomial(polynomial)
        converted = OccupationPolynomial{D,S}([
            monomial => convert(D, value) for (monomial, value) in occupation
        ])
        _push_finite_width_occupation!(out, term, converted)
    end
    return out
end

"""
    finite_width_occupation_collision(collision)

Lower the analytically resolved branch of a `FiniteWidthEvaluatedCollision` with the same exact
statistics-generic rules as the strict-quasiparticle collision compiler,

```math
F_S = 1 + 2σ_S n_S,
C_n = (σ_S / 2) I_package.
```

The affine Kadanoff--Baym pieces are recombined at this stage: the external target distribution
factor is included only for the original distribution-coefficient branch. Unsupported finite-width
frequency structures remain explicit as separate offset/distribution maps and are not projected to
occupations.

The Lorentzian spectral model used to obtain `collision` is treated as fixed. A self-consistent
linear response must additionally differentiate any state-dependent linewidth `Γ[n]`; that chain
rule is deliberately outside this operation.
"""
function finite_width_occupation_collision(
    collision::FiniteWidthEvaluatedCollision{C,S,O,E1,E2,G,Ctx}
) where {C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    D = promote_type(C, Rational{Int})
    terms = Dict{FiniteWidthFrequencyTerm{S,E1,E2},OccupationPolynomial{D,S}}()
    target = target_family(collision)

    _append_finite_width_occupation_component!(
        terms, collision_offset(collision), target, false
    )
    _append_finite_width_occupation_component!(
        terms, collision_distribution_coefficient(collision), target, true
    )

    unsupported_offset = Dict{FiniteWidthFrequencyTerm{S,E1,E2},D}(
        term => convert(D, coefficient) for
        (term, coefficient) in finite_width_unsupported_terms(collision_offset(collision))
    )
    unsupported_distribution = Dict{FiniteWidthFrequencyTerm{S,E1,E2},D}(
        term => convert(D, coefficient) for (term, coefficient) in
        finite_width_unsupported_terms(collision_distribution_coefficient(collision))
    )

    return FiniteWidthOccupationCollision{D,S,O,E1,E2,G,Ctx}(
        terms,
        unsupported_offset,
        unsupported_distribution,
        target,
        parameters(collision),
        wigner_context(collision),
    )
end
