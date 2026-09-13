@enum CausalFrequencyIntegrationKind::UInt8 begin
    CausalFrequencyIntegrated = 0
    CausalFrequencyMarginalAtInfinity = 1
    CausalFrequencyDivergentAtInfinity = 2
    CausalFrequencyZeroPrescription = 3
    CausalFrequencyRepeatedPole = 4
end

"""
One deterministic expression-level causal-frequency integration step.

`expression` is either the integrated canonical expression or, for a deferred result, the
unchanged canonical input promoted to the coefficient type used by causal residues. `kind`
records why a frequency could not yet be integrated. Deferred cases are semantic outputs for
later reduction layers, not internal errors.
"""
struct CausalFrequencyIntegrationResult{C<:Number,S<:Statistics}
    expression::CausalFrequencyExpression{C,S}
    kind::CausalFrequencyIntegrationKind
    frequency_index::Int
end

causal_frequency_integration_kind(result::CausalFrequencyIntegrationResult) = result.kind
function causal_frequency_integration_expression(result::CausalFrequencyIntegrationResult)
    return result.expression
end
function causal_frequency_integration_index(result::CausalFrequencyIntegrationResult)
    return result.frequency_index
end

"""
A statically sized order in which loop frequencies are reduced.

The plan stores only distinct positive frequency-basis indices. It does not encode any physics
choice: different complete plans for the same regular expression must give the same canonical
result. The canonical plan is ascending basis order; alternative plans exist to certify this
order-invariance contract explicitly.
"""
struct CausalFrequencyReductionPlan{N}
    order::NTuple{N,Int}

    function CausalFrequencyReductionPlan(order::NTuple{N,Int}) where {N}
        all(index -> index > 0, order) ||
            throw(ArgumentError("causal frequency-reduction indices must be positive"))
        length(unique(order)) == N ||
            throw(ArgumentError("causal frequency-reduction indices must be distinct"))
        return new{N}(order)
    end
end

causal_frequency_reduction_order(plan::CausalFrequencyReductionPlan) = plan.order
function canonical_causal_frequency_reduction_plan(::Val{N}) where {N}
    return CausalFrequencyReductionPlan(ntuple(identity, Val(N)))
end

"""
Result of applying a complete or partial causal-frequency reduction plan.

`completed_steps` counts successful contour integrations. If `kind` is deferred,
`blocked_frequency` is the first basis index whose value belongs to a later semantic layer.
For a completed plan it is zero.
"""
struct CausalFrequencyReductionResult{C<:Number,S<:Statistics,N}
    expression::CausalFrequencyExpression{C,S}
    plan::CausalFrequencyReductionPlan{N}
    completed_steps::Int
    kind::CausalFrequencyIntegrationKind
    blocked_frequency::Int
end

causal_frequency_reduction_expression(result::CausalFrequencyReductionResult) = result.expression
causal_frequency_reduction_plan(result::CausalFrequencyReductionResult) = result.plan
causal_frequency_reduction_steps(result::CausalFrequencyReductionResult) = result.completed_steps
causal_frequency_reduction_kind(result::CausalFrequencyReductionResult) = result.kind
function causal_frequency_reduction_blocked_frequency(result::CausalFrequencyReductionResult)
    return result.blocked_frequency
end

function _complex_causal_frequency_expression(
    expression::CausalFrequencyExpression{C,S}
) where {C<:Number,S<:Statistics}
    D = promote_type(C, ComplexRationals)
    terms = CausalFrequencyTerm{D,S}[]
    sizehint!(terms, length(expression.terms))
    for term in expression.terms
        push!(terms, CausalFrequencyTerm(convert(D, term.coefficient), term.denominators))
    end
    return CausalFrequencyExpression(terms)
end

function _substitute_causal_frequency_pole(
    denominator::CausalFrequencyDenominator{S},
    pivot::CausalFrequencyDenominator{S},
    frequency_index::Int,
) where {S<:Statistics}
    pivot_coefficient = pivot.loop_coefficients[frequency_index]
    iszero(pivot_coefficient) &&
        throw(ArgumentError("causal residue pivot does not depend on selected frequency"))
    coefficient = denominator.loop_coefficients[frequency_index]
    iszero(coefficient) && return denominator

    ratio = coefficient / pivot_coefficient
    transformed = copy(denominator.loop_coefficients)
    @inbounds for i in eachindex(transformed, pivot.loop_coefficients)
        transformed[i] -= ratio * pivot.loop_coefficients[i]
    end
    transformed[frequency_index] = zero(EnergyCoefficient)
    return CausalFrequencyDenominator{S}(
        transformed,
        denominator.energy - ratio * pivot.energy,
        denominator.infinitesimal - ratio * pivot.infinitesimal,
    )
end

@inline function _causal_pole_is_upper(
    denominator::CausalFrequencyDenominator, frequency_index::Int
)
    coefficient = denominator.loop_coefficients[frequency_index]
    iszero(coefficient) &&
        throw(ArgumentError("causal pole query requires a frequency-dependent denominator"))
    iszero(denominator.infinitesimal) &&
        throw(ArgumentError("zero prescription has no causal half-plane"))
    return denominator.infinitesimal / coefficient < 0
end

@inline function _identically_zero_causal_denominator(
    denominator::CausalFrequencyDenominator
)
    return all(iszero, denominator.loop_coefficients) &&
           iszero(denominator.energy) &&
           iszero(denominator.infinitesimal)
end

function _causal_frequency_blocker(
    expression::CausalFrequencyExpression, frequency_index::Int
)::CausalFrequencyIntegrationKind
    for term in expression.terms
        denominators = term.denominators
        for pivot_index in eachindex(denominators)
            pivot = denominators[pivot_index]
            iszero(pivot.loop_coefficients[frequency_index]) && continue
            iszero(pivot.infinitesimal) && return CausalFrequencyZeroPrescription
            _causal_pole_is_upper(pivot, frequency_index) || continue

            for denominator_index in eachindex(denominators)
                denominator_index == pivot_index && continue
                transformed = _substitute_causal_frequency_pole(
                    denominators[denominator_index], pivot, frequency_index
                )
                _identically_zero_causal_denominator(transformed) &&
                    return CausalFrequencyRepeatedPole
            end
        end
    end
    return CausalFrequencyIntegrated
end

function _causal_frequency_residues(
    term::CausalFrequencyTerm{C,S}, frequency_index::Int
) where {C<:Number,S<:Statistics}
    D = promote_type(C, ComplexRationals)
    out = CausalFrequencyTerm{D,S}[]
    for pivot_index in eachindex(term.denominators)
        pivot = term.denominators[pivot_index]
        pivot_coefficient = pivot.loop_coefficients[frequency_index]
        iszero(pivot_coefficient) && continue
        _causal_pole_is_upper(pivot, frequency_index) || continue

        denominators = CausalFrequencyDenominator{S}[]
        sizehint!(denominators, length(term.denominators) - 1)
        for denominator_index in eachindex(term.denominators)
            denominator_index == pivot_index && continue
            push!(
                denominators,
                _substitute_causal_frequency_pole(
                    term.denominators[denominator_index], pivot, frequency_index
                ),
            )
        end

        residue_factor = convert(D, im) / convert(D, pivot_coefficient)
        push!(
            out,
            CausalFrequencyTerm(
                convert(D, term.coefficient) * residue_factor, denominators
            ),
        )
    end
    return out
end

"""
    integrate_causal_frequency_expression(expression, frequency_index)

Integrate one loop frequency only after the *complete canonical expression* is certified to
have a vanishing large semicircle. The contour is closed deterministically in the upper half
plane with measure `dω/(2π)`.

This differs essentially from termwise contour reduction: individual summands may have a
marginal `1/ω` tail. Their residues are still linear contributions to the integral whenever the
sum has `O(1/ω²)` decay. The global decay test therefore precedes every per-term residue
calculation.

Geometries whose value belongs to a later semantic layer are returned unchanged with a typed
status: marginal or divergent behavior at infinity, a frequency-dependent zero prescription,
or an upper-half-plane repeated pole. No physical value is assigned to those cases here.
"""
function integrate_causal_frequency_expression(
    expression::CausalFrequencyExpression{C,S}, frequency_index::Int
) where {C<:Number,S<:Statistics}
    canonical = _complex_causal_frequency_expression(expression)
    isempty(canonical) && return CausalFrequencyIntegrationResult(
        canonical, CausalFrequencyIntegrated, frequency_index
    )

    # A nonzero denominator-free term is a genuine O(1) contribution. The canonical
    # constructor has already combined all constant terms, so reaching this branch proves that
    # the complete expression does not decay in any remaining loop frequency. Classify it
    # before denominator-basis validation: UV/infinity behavior is a semantic result, not an
    # invalid-expression exception.
    if any(term -> isempty(term.denominators), canonical.terms)
        return CausalFrequencyIntegrationResult(
            canonical, CausalFrequencyDivergentAtInfinity, frequency_index
        )
    end

    _check_causal_frequency_index(canonical, frequency_index)
    decay = frequency_decay_lower_bound(canonical, frequency_index)
    if decay < 1
        return CausalFrequencyIntegrationResult(
            canonical, CausalFrequencyDivergentAtInfinity, frequency_index
        )
    elseif decay == 1
        return CausalFrequencyIntegrationResult(
            canonical, CausalFrequencyMarginalAtInfinity, frequency_index
        )
    end

    blocker = _causal_frequency_blocker(canonical, frequency_index)
    blocker === CausalFrequencyIntegrated ||
        return CausalFrequencyIntegrationResult(canonical, blocker, frequency_index)

    D = promote_type(C, ComplexRationals)
    residues = CausalFrequencyTerm{D,S}[]
    for term in canonical.terms
        append!(residues, _causal_frequency_residues(term, frequency_index))
    end
    integrated = CausalFrequencyExpression(residues)
    return CausalFrequencyIntegrationResult(
        integrated, CausalFrequencyIntegrated, frequency_index
    )
end

"""
    reduce_causal_frequency_expression(expression, plan)

Apply expression-level contour integration in the explicit order encoded by `plan`, stopping at
the first typed deferred geometry. The canonical plan is ascending loop-frequency basis order,
but no physical result is allowed to depend on that choice: regular expressions must reduce to
the same canonical result under any complete admissible plan.
"""
function reduce_causal_frequency_expression(
    expression::CausalFrequencyExpression{C,S}, plan::CausalFrequencyReductionPlan{N}
) where {C<:Number,S<:Statistics,N}
    current = _complex_causal_frequency_expression(expression)
    completed_steps = 0
    for frequency_index in plan.order
        step = integrate_causal_frequency_expression(current, frequency_index)
        current = causal_frequency_integration_expression(step)
        kind = causal_frequency_integration_kind(step)
        kind === CausalFrequencyIntegrated ||
            return CausalFrequencyReductionResult(
                current, plan, completed_steps, kind, frequency_index
            )
        completed_steps += 1
    end
    return CausalFrequencyReductionResult(
        current, plan, completed_steps, CausalFrequencyIntegrated, 0
    )
end
