"""One causal contribution together with its pre-pivot real affine support provenance."""
struct CausalFrequencySupportContribution{C<:Number,S<:Statistics}
    term::CausalFrequencyTerm{C,S}
    support::AffineSingularSupport{S}
end

function _support_contribution_isless(
    a::CausalFrequencySupportContribution{C,S}, b::CausalFrequencySupportContribution{C,S}
) where {C<:Number,S<:Statistics}
    return _causal_denominator_vectors_isless(a.term.denominators, b.term.denominators)
end

function _causal_expression_from_support_contributions(
    contributions::Vector{CausalFrequencySupportContribution{C,S}}
) where {C<:Number,S<:Statistics}
    terms = CausalFrequencyTerm{C,S}[contribution.term for contribution in contributions]
    return CausalFrequencyExpression(terms)
end

function _prune_cancelled_support_contributions(
    contributions::Vector{CausalFrequencySupportContribution{C,S}}
) where {C<:Number,S<:Statistics}
    isempty(contributions) && return contributions
    ordered = copy(contributions)
    sort!(ordered; lt=_support_contribution_isless)

    out = CausalFrequencySupportContribution{C,S}[]
    sizehint!(out, length(ordered))
    first_index = firstindex(ordered)
    while first_index <= lastindex(ordered)
        denominators = ordered[first_index].term.denominators
        last_index = first_index
        coefficient = zero(C)
        while last_index <= lastindex(ordered) &&
              isequal(ordered[last_index].term.denominators, denominators)
            coefficient += ordered[last_index].term.coefficient
            last_index += 1
        end
        if !iszero(coefficient)
            append!(out, @view ordered[first_index:(last_index - 1)])
        end
        first_index = last_index
    end
    return out
end

function _initial_support_contributions(
    expression::CausalFrequencyExpression{C,S}
) where {C<:Number,S<:Statistics}
    canonical = _complex_causal_frequency_expression(expression)
    D = promote_type(C, ComplexRationals)
    contributions = CausalFrequencySupportContribution{D,S}[]
    sizehint!(contributions, length(canonical.terms))
    for term in canonical.terms
        push!(
            contributions,
            CausalFrequencySupportContribution{D,S}(term, affine_singular_support(term)),
        )
    end
    return contributions
end

function _propagate_support_residues(
    contributions::Vector{CausalFrequencySupportContribution{C,S}}, frequency_index::Int
) where {C<:Number,S<:Statistics}
    out = CausalFrequencySupportContribution{C,S}[]
    for contribution in contributions
        for residue in _causal_frequency_residues(contribution.term, frequency_index)
            push!(
                out, CausalFrequencySupportContribution{C,S}(residue, contribution.support)
            )
        end
    end
    return _prune_cancelled_support_contributions(out)
end

"""
Result of causal-frequency reduction with exact pre-pivot affine-support provenance.

The support labels never participate in contour-safety decisions. `expression` is the complete
canonical sum reconstructed from `contributions`; frozen #299 remains the sole authority for
infinity, pole and blocker semantics.
"""
struct CausalFrequencySupportReductionResult{C<:Number,S<:Statistics,N}
    expression::CausalFrequencyExpression{C,S}
    contributions::Vector{CausalFrequencySupportContribution{C,S}}
    plan::CausalFrequencyReductionPlan{N}
    completed_steps::Int
    kind::CausalFrequencyIntegrationKind
    blocked_frequency::Int
end

function causal_frequency_support_expression(result::CausalFrequencySupportReductionResult)
    return result.expression
end
function causal_frequency_support_kind(result::CausalFrequencySupportReductionResult)
    return result.kind
end
function causal_frequency_support_steps(result::CausalFrequencySupportReductionResult)
    return result.completed_steps
end
function causal_frequency_support_blocked_frequency(
    result::CausalFrequencySupportReductionResult
)
    return result.blocked_frequency
end

"""
    reduce_causal_frequency_with_support(expression, plan)

Apply the frozen #299 contour semantics while carrying affine singular-support provenance beside
each linear contribution. Before every integration step, all contributions are recombined into
the complete canonical expression and that complete sum is passed to
`integrate_causal_frequency_expression`. Only a globally certified regular step is then
propagated contribution-by-contribution through the same exact residue kernel.

This preserves cancellations of marginal tails and residue contributions across different
support histories. If one canonical denominator product cancels exactly in the complete sum,
all provenance attached solely to that zero contribution is discarded.
"""
function reduce_causal_frequency_with_support(
    expression::CausalFrequencyExpression{C,S}, plan::CausalFrequencyReductionPlan{N}
) where {C<:Number,S<:Statistics,N}
    contributions = _initial_support_contributions(expression)
    current = _causal_expression_from_support_contributions(contributions)
    completed_steps = 0

    for frequency_index in plan.order
        step = integrate_causal_frequency_expression(current, frequency_index)
        kind = causal_frequency_integration_kind(step)
        if kind !== CausalFrequencyIntegrated
            return CausalFrequencySupportReductionResult(
                current, contributions, plan, completed_steps, kind, frequency_index
            )
        end

        contributions = _propagate_support_residues(contributions, frequency_index)
        current = _causal_expression_from_support_contributions(contributions)
        expected = causal_frequency_integration_expression(step)
        current == expected ||
            throw(ErrorException("support provenance changed canonical causal reduction"))
        completed_steps += 1
    end

    return CausalFrequencySupportReductionResult(
        current, contributions, plan, completed_steps, CausalFrequencyIntegrated, 0
    )
end

"""
Return the distinct pre-pivot affine supports whose current contributions realize the typed
singular blocker. Boundary-value prescriptions remain in the causal terms and are not folded
into physical support identity here.
"""
function causal_frequency_singular_supports(
    result::CausalFrequencySupportReductionResult{C,S}
) where {C<:Number,S<:Statistics}
    kind = result.kind
    (kind === CausalFrequencyPinch || kind === CausalFrequencyRepeatedPole) ||
        return AffineSingularSupport{S}[]

    supports = AffineSingularSupport{S}[]
    for contribution in result.contributions
        local_expression = CausalFrequencyExpression(
            CausalFrequencyTerm{C,S}[contribution.term]
        )
        _causal_frequency_blocker(local_expression, result.blocked_frequency) === kind ||
            continue
        support = contribution.support
        has_affine_singular_support(support) || continue
        any(existing -> isequal(existing, support), supports) || push!(supports, support)
    end
    return supports
end
