"""
One affine finite-width collision component after explicit Lorentzian spectral-model evaluation.

`resolved` stores the original classified finite-width term together with its coefficient after
multiplying by the exact factorized or residual-convolution spectral weight. `unsupported` keeps
terms for which the analytic finite-width backend has no value; their original collision
coefficient is retained without an implicit cutoff or sharp-shell fallback.
"""
struct FiniteWidthEvaluatedExpression{
    C<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext
}
    resolved::Dict{FiniteWidthFrequencyTerm{S,E1,E2},C}
    unsupported::Dict{FiniteWidthFrequencyTerm{S,E1,E2},C}
    context::Ctx
end

statistics(::FiniteWidthEvaluatedExpression{C,S}) where {C,S} = S
gradient_order(::FiniteWidthEvaluatedExpression{C,S,E1,E2,G}) where {C,S,E1,E2,G} = Val(G)
wigner_context(expression::FiniteWidthEvaluatedExpression) = expression.context
function finite_width_resolved_terms(expression::FiniteWidthEvaluatedExpression)
    return expression.resolved
end
function finite_width_unsupported_terms(expression::FiniteWidthEvaluatedExpression)
    return expression.unsupported
end
function Base.isempty(expression::FiniteWidthEvaluatedExpression)
    return isempty(expression.resolved) && isempty(expression.unsupported)
end

function _evaluate_finite_width_expression(
    expression::FiniteWidthFrequencyExpression{C,S,E1,E2,G,Ctx},
    model::LorentzianSpectralModel{M,S},
    ω_external::A,
) where {C<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext,M<:Real,A<:Real}
    D = promote_type(C, M, A)
    resolved = Dict{FiniteWidthFrequencyTerm{S,E1,E2},D}()
    unsupported = Dict{FiniteWidthFrequencyTerm{S,E1,E2},D}()
    for (term, coefficient) in expression
        kind = finite_width_reduction_kind(term)
        if kind === FiniteWidthFactorized
            weight = evaluate_spectral_weight(finite_width_spectral_reduction(term), model)
            resolved[term] = convert(D, coefficient) * convert(D, weight)
        elseif kind === FiniteWidthConvolution
            weight = evaluate_spectral_convolution(
                finite_width_convolution_reduction(term), model, ω_external
            )
            resolved[term] = convert(D, coefficient) * convert(D, weight)
        else
            unsupported[term] = convert(D, coefficient)
        end
    end
    return FiniteWidthEvaluatedExpression{D,S,E1,E2,G,Ctx}(
        resolved, unsupported, wigner_context(expression)
    )
end

"""
Complete affine collision after evaluating all analytically supported finite-width frequency
weights against an explicit Lorentzian spectral model.

The two Kadanoff--Baym components remain separate. This type still carries the original
statistical-weight tags through its `FiniteWidthFrequencyTerm` keys; it does not perform
`F -> n`, momentum integration, a trap closure, or collective-mode projection.
"""
struct FiniteWidthEvaluatedCollision{
    C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext
}
    offset::FiniteWidthEvaluatedExpression{C,S,E1,E2,G,Ctx}
    distribution_coefficient::FiniteWidthEvaluatedExpression{C,S,E1,E2,G,Ctx}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::FiniteWidthEvaluatedCollision{C,S,O}) where {C,S,O} = O
statistics(::FiniteWidthEvaluatedCollision{C,S}) where {C,S} = S
parameters(collision::FiniteWidthEvaluatedCollision) = collision.parameter
function gradient_order(
    ::FiniteWidthEvaluatedCollision{C,S,O,E1,E2,G}
) where {C,S,O,E1,E2,G}
    return Val(G)
end
wigner_context(collision::FiniteWidthEvaluatedCollision) = collision.context
target_family(collision::FiniteWidthEvaluatedCollision) = collision.target
collision_offset(collision::FiniteWidthEvaluatedCollision) = collision.offset
function collision_distribution_coefficient(collision::FiniteWidthEvaluatedCollision)
    return collision.distribution_coefficient
end

"""
    evaluate_finite_width_collision(collision, model, ω_external)

Evaluate every analytically resolved term of a `FiniteWidthFrequencyCollision` using the supplied
`LorentzianSpectralModel`. Factorized terms use the exact integrated Lorentzian powers, while
one-residual-shell terms use the exact Cauchy convolution at the explicitly supplied microscopic
spectral frequency `ω_external`.

Unsupported terms remain explicit and unevaluated. Missing spectral data are errors from the
underlying spectral model; this operation never inserts a cutoff, zero width, or strict-QP value.
The external spectral frequency here is not the center-time collective response frequency.
"""
function evaluate_finite_width_collision(
    collision::FiniteWidthFrequencyCollision{C,S,O,E1,E2,G,Ctx},
    model::LorentzianSpectralModel{M,S},
    ω_external::A,
) where {C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext,M<:Real,A<:Real}
    offset = _evaluate_finite_width_expression(
        collision_offset(collision), model, ω_external
    )
    distribution = _evaluate_finite_width_expression(
        collision_distribution_coefficient(collision), model, ω_external
    )
    D = promote_type(C, M, A)
    return FiniteWidthEvaluatedCollision{D,S,O,E1,E2,G,Ctx}(
        offset,
        distribution,
        target_family(collision),
        parameters(collision),
        wigner_context(collision),
    )
end
