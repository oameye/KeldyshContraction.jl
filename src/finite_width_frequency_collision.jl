"""Structural finite-width frequency reduction selected for one collision term."""
@enum FiniteWidthFrequencyReductionKind::UInt8 begin
    FiniteWidthFactorized = 0
    FiniteWidthConvolution = 1
    FiniteWidthUnsupported = 2
end

"""
One spectral/dispersive collision term together with its exact finite-width frequency reduction.

The original `SpectralDispersiveTerm` is retained verbatim so statistical-weight tags, exact
momentum routing, derivative kinematics, equal-time shifts, and graph provenance remain
inspectable. The external target family is part of the identity because the residual-convolution
mismatch depends on it. `FiniteWidthFactorized` stores the repeated/independent-line Lorentzian
reduction, `FiniteWidthConvolution` stores the one-residual-shell Cauchy reduction, and
`FiniteWidthUnsupported` keeps the source term without assigning a finite value.
"""
struct FiniteWidthFrequencyTerm{S<:Statistics,E1,E2}
    source::SpectralDispersiveTerm{S,E1,E2}
    target::FieldFamily{S}
    kind::FiniteWidthFrequencyReductionKind
    spectral::LorentzianSpectralReduction{S}
    convolution::LorentzianConvolutionReduction{S}
end

statistics(::FiniteWidthFrequencyTerm{S}) where {S<:Statistics} = S
target_family(term::FiniteWidthFrequencyTerm) = term.target
finite_width_source_term(term::FiniteWidthFrequencyTerm) = term.source
finite_width_reduction_kind(term::FiniteWidthFrequencyTerm) = term.kind

function finite_width_spectral_reduction(term::FiniteWidthFrequencyTerm)
    term.kind === FiniteWidthFactorized ||
        throw(ArgumentError("finite-width term is not a factorized Lorentzian reduction"))
    return term.spectral
end

function finite_width_convolution_reduction(term::FiniteWidthFrequencyTerm)
    term.kind === FiniteWidthConvolution ||
        throw(ArgumentError("finite-width term is not a residual Lorentzian convolution"))
    return term.convolution
end

function Base.isequal(
    a::FiniteWidthFrequencyTerm{S,E1,E2}, b::FiniteWidthFrequencyTerm{S,E1,E2}
) where {S,E1,E2}
    return a.kind === b.kind && isequal(a.target, b.target) && isequal(a.source, b.source)
end
Base.:(==)(a::FiniteWidthFrequencyTerm, b::FiniteWidthFrequencyTerm) = isequal(a, b)
function Base.hash(term::FiniteWidthFrequencyTerm, h::UInt)
    return hash(
        term.kind, hash(term.target, hash(term.source, hash(FiniteWidthFrequencyTerm, h)))
    )
end

function _finite_width_frequency_term(
    term::SpectralDispersiveTerm{S,E1,E2}, target::FieldFamily{S}
) where {S<:Statistics,E1,E2}
    spectral = lorentzian_spectral_reduction(term)
    convolution = lorentzian_convolution_reduction(term, target)
    kind = if spectral_reduction_resolved(spectral)
        FiniteWidthFactorized
    elseif convolution_reduction_resolved(convolution)
        FiniteWidthConvolution
    else
        FiniteWidthUnsupported
    end
    return FiniteWidthFrequencyTerm{S,E1,E2}(term, target, kind, spectral, convolution)
end

"""Exact finite-width frequency classification of one affine collision component."""
struct FiniteWidthFrequencyExpression{
    C<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext
}
    terms::Dict{FiniteWidthFrequencyTerm{S,E1,E2},C}
    context::Ctx
end

statistics(::FiniteWidthFrequencyExpression{C,S}) where {C,S} = S
gradient_order(::FiniteWidthFrequencyExpression{C,S,E1,E2,G}) where {C,S,E1,E2,G} = Val(G)
wigner_context(expression::FiniteWidthFrequencyExpression) = expression.context
Base.length(expression::FiniteWidthFrequencyExpression) = length(expression.terms)
Base.isempty(expression::FiniteWidthFrequencyExpression) = isempty(expression.terms)
Base.iszero(expression::FiniteWidthFrequencyExpression) = isempty(expression.terms)
Base.iterate(expression::FiniteWidthFrequencyExpression) = iterate(expression.terms)
function Base.iterate(expression::FiniteWidthFrequencyExpression, state)
    return iterate(expression.terms, state)
end
function Base.eltype(
    ::Type{FiniteWidthFrequencyExpression{C,S,E1,E2,G,Ctx}}
) where {C,S,E1,E2,G,Ctx}
    return Pair{FiniteWidthFrequencyTerm{S,E1,E2},C}
end

function _finite_width_frequency_expression(
    expression::SpectralDispersiveExpression{C,S,E1,E2,G,Ctx}, target::FieldFamily{S}
) where {C<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext}
    terms = Dict{FiniteWidthFrequencyTerm{S,E1,E2},C}()
    sizehint!(terms, length(expression))
    for (source, coefficient) in expression
        terms[_finite_width_frequency_term(source, target)] = coefficient
    end
    return FiniteWidthFrequencyExpression{C,S,E1,E2,G,Ctx}(
        terms, wigner_context(expression)
    )
end

"""
Collision-level finite-width frequency classification before spectral-model evaluation.

The exact affine Kadanoff--Baym structure

    I_coll = I_0 + F_target(k) I_1

is retained as two `FiniteWidthFrequencyExpression`s. This representation performs no spectral
model evaluation, no external-frequency choice, no occupation substitution, and no momentum
integration. Unsupported/dispersive/Trotter structures remain attached to their original source
terms through `FiniteWidthUnsupported`.
"""
struct FiniteWidthFrequencyCollision{
    C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext
}
    offset::FiniteWidthFrequencyExpression{C,S,E1,E2,G,Ctx}
    distribution_coefficient::FiniteWidthFrequencyExpression{C,S,E1,E2,G,Ctx}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::FiniteWidthFrequencyCollision{C,S,O}) where {C,S,O} = O
statistics(::FiniteWidthFrequencyCollision{C,S}) where {C,S} = S
parameters(collision::FiniteWidthFrequencyCollision) = collision.parameter
function gradient_order(
    ::FiniteWidthFrequencyCollision{C,S,O,E1,E2,G}
) where {C,S,O,E1,E2,G}
    return Val(G)
end
wigner_context(collision::FiniteWidthFrequencyCollision) = collision.context
target_family(collision::FiniteWidthFrequencyCollision) = collision.target
collision_offset(collision::FiniteWidthFrequencyCollision) = collision.offset
function collision_distribution_coefficient(collision::FiniteWidthFrequencyCollision)
    return collision.distribution_coefficient
end

"""
    finite_width_frequency_collision(collision)

Classify every term of a complete `SpectralDispersiveCollision` with the existing exact
Lorentzian frequency backends. Factorized repeated-line sectors and one-residual-shell
convolutions become explicit resolved terms; all other structures remain explicit unsupported
terms. The original source term and affine collision component are preserved exactly.

This stage only classifies frequency structure. In particular it does not supply linewidth data,
evaluate a microscopic external spectral frequency, lower `F -> n`, or perform loop-momentum
integration.
"""
function finite_width_frequency_collision(
    collision::SpectralDispersiveCollision{C,S,O,E1,E2,G,Ctx}
) where {C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    return FiniteWidthFrequencyCollision{C,S,O,E1,E2,G,Ctx}(
        _finite_width_frequency_expression(
            collision_offset(collision), target_family(collision)
        ),
        _finite_width_frequency_expression(
            collision_distribution_coefficient(collision), target_family(collision)
        ),
        target_family(collision),
        parameters(collision),
        wigner_context(collision),
    )
end
