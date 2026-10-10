"""
Explicit occupation-channel response of a Lorentzian spectral model.

Each entry maps one kinetic variation channel `δn_a` to the line-resolved first variation
`(δE, δΓ, δZ)` of the supplied spectral model. An omitted occupation channel means that the
spectral model has no response in that channel. For every channel that is supplied, the
`LorentzianSpectralModelVariation` must still contain explicit variations, including explicit
zeros, for every spectral line required by a resolved finite-width term.
"""
struct OccupationSpectralResponse{V<:Real,S<:Statistics}
    terms::Vector{Pair{OccupationAtom{S},LorentzianSpectralModelVariation{V,S}}}

    function OccupationSpectralResponse{V,S}(
        terms::Vector{Pair{OccupationAtom{S},LorentzianSpectralModelVariation{V,S}}},
        ::Val{:raw},
    ) where {V<:Real,S<:Statistics}
        return new{V,S}(terms)
    end
end

function OccupationSpectralResponse{V,S}(
    terms::Vector{Pair{OccupationAtom{S},LorentzianSpectralModelVariation{V,S}}}
) where {V<:Real,S<:Statistics}
    ordered = sort(copy(terms); by=first)
    for index in 2:length(ordered)
        isequal(first(ordered[index - 1]), first(ordered[index])) && throw(
            ArgumentError(
                "spectral occupation response contains a duplicate variation channel"
            ),
        )
    end
    return OccupationSpectralResponse{V,S}(ordered, Val(:raw))
end

function OccupationSpectralResponse(
    terms::AbstractDict{OccupationAtom{S},LorentzianSpectralModelVariation{V,S}}
) where {V<:Real,S<:Statistics}
    return OccupationSpectralResponse{V,S}(collect(terms))
end

Base.length(response::OccupationSpectralResponse) = length(response.terms)
Base.isempty(response::OccupationSpectralResponse) = isempty(response.terms)
Base.iterate(response::OccupationSpectralResponse) = iterate(response.terms)
function Base.iterate(response::OccupationSpectralResponse, state)
    return iterate(response.terms, state)
end
function Base.eltype(::Type{OccupationSpectralResponse{V,S}}) where {V,S}
    return Pair{OccupationAtom{S},LorentzianSpectralModelVariation{V,S}}
end
Base.IteratorSize(::Type{<:OccupationSpectralResponse}) = Base.HasLength()

"""Return the explicit occupation-channel to spectral-model-variation mapping."""
spectral_occupation_response_terms(response::OccupationSpectralResponse) = response.terms

function _finite_width_weight(
    term::FiniteWidthFrequencyTerm{S}, model::LorentzianSpectralModel{M,S}, ω_external::A
) where {S<:Statistics,M<:Real,A<:Real}
    kind = finite_width_reduction_kind(term)
    if kind === FiniteWidthFactorized
        return evaluate_spectral_weight(finite_width_spectral_reduction(term), model)
    elseif kind === FiniteWidthConvolution
        return evaluate_spectral_convolution(
            finite_width_convolution_reduction(term), model, ω_external
        )
    end
    return throw(ArgumentError("unsupported finite-width term has no spectral weight"))
end

function _finite_width_weight_variation(
    term::FiniteWidthFrequencyTerm{S},
    model::LorentzianSpectralModel{M,S},
    variation::LorentzianSpectralModelVariation{V,S},
    ω_external::A,
) where {S<:Statistics,M<:Real,V<:Real,A<:Real}
    kind = finite_width_reduction_kind(term)
    if kind === FiniteWidthFactorized
        return evaluate_spectral_weight_variation(
            finite_width_spectral_reduction(term), model, variation
        )
    elseif kind === FiniteWidthConvolution
        return evaluate_spectral_convolution_variation(
            finite_width_convolution_reduction(term), model, variation, ω_external
        )
    end
    return throw(
        ArgumentError("unsupported finite-width term has no spectral-weight variation")
    )
end

function _finite_width_occupation_polynomial(
    term::FiniteWidthFrequencyTerm{S},
    coefficient::C,
    target::FieldFamily{S},
    include_external::Bool,
) where {C<:Number,S<:Statistics}
    statistical = _collision_statistical_monomial(
        finite_width_source_term(term), target, include_external
    )
    polynomial = StatisticalPolynomial{C,S}([statistical => coefficient])
    return occupation_collision_polynomial(polynomial)
end

function _scale_background_linearization(
    linearization::BackgroundOccupationLinearization{C,S}, scale::A, ::Type{D}
) where {C<:Number,A<:Number,D<:Number,S<:Statistics}
    factor = convert(D, scale)
    return BackgroundOccupationLinearization{D,S}([
        atom => factor * convert(D, coefficient) for (atom, coefficient) in linearization
    ])
end

function _combine_background_linearizations(
    first_linearization::BackgroundOccupationLinearization{A,S},
    second_linearization::BackgroundOccupationLinearization{B,S},
    ::Type{D},
) where {A<:Number,B<:Number,D<:Number,S<:Statistics}
    terms = Pair{OccupationAtom{S},D}[]
    sizehint!(terms, length(first_linearization) + length(second_linearization))
    for (atom, coefficient) in first_linearization
        push!(terms, atom => convert(D, coefficient))
    end
    for (atom, coefficient) in second_linearization
        push!(terms, atom => convert(D, coefficient))
    end
    return BackgroundOccupationLinearization{D,S}(terms)
end

function _push_background_linearization!(
    out::Dict{FiniteWidthFrequencyTerm{S,E1,E2},BackgroundOccupationLinearization{D,S}},
    term::FiniteWidthFrequencyTerm{S,E1,E2},
    linearization::BackgroundOccupationLinearization{D,S},
    ::Type{D},
) where {D<:Number,S<:Statistics,E1,E2}
    iszero(linearization) && return out
    if haskey(out, term)
        combined = _combine_background_linearizations(out[term], linearization, D)
        if iszero(combined)
            delete!(out, term)
        else
            out[term] = combined
        end
    else
        out[term] = linearization
    end
    return out
end

"""
Self-consistent finite-width collision response around a supplied occupation background.

For every analytically resolved finite-width term, the response is retained in three inspectable
forms: the fixed-model occupation response, the spectral-data response, and their exact product-rule
sum.

`occupation` stores the first contribution, `spectral` stores the second, and `terms` stores their
exact sum. All three keep the original `FiniteWidthFrequencyTerm` key. Unsupported affine sectors
remain explicit and are not assigned a response value.
"""
struct FiniteWidthCollisionLinearization{
    C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext
}
    terms::Dict{FiniteWidthFrequencyTerm{S,E1,E2},BackgroundOccupationLinearization{C,S}}
    occupation::Dict{
        FiniteWidthFrequencyTerm{S,E1,E2},BackgroundOccupationLinearization{C,S}
    }
    spectral::Dict{FiniteWidthFrequencyTerm{S,E1,E2},BackgroundOccupationLinearization{C,S}}
    unsupported_offset::Dict{FiniteWidthFrequencyTerm{S,E1,E2},C}
    unsupported_distribution::Dict{FiniteWidthFrequencyTerm{S,E1,E2},C}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::FiniteWidthCollisionLinearization{C,S,O}) where {C,S,O} = O
statistics(::FiniteWidthCollisionLinearization{C,S}) where {C,S} = S
parameters(response::FiniteWidthCollisionLinearization) = response.parameter
target_family(response::FiniteWidthCollisionLinearization) = response.target
function gradient_order(
    ::FiniteWidthCollisionLinearization{C,S,O,E1,E2,G}
) where {C,S,O,E1,E2,G}
    return Val(G)
end
wigner_context(response::FiniteWidthCollisionLinearization) = response.context

"""Return the complete resolved finite-width collision linear-response terms."""
finite_width_linearized_terms(response::FiniteWidthCollisionLinearization) = response.terms

"""Return the fixed-spectral occupation-response contribution terms."""
function finite_width_occupation_response_terms(response::FiniteWidthCollisionLinearization)
    return response.occupation
end

"""Return the spectral-data response contribution terms."""
function finite_width_spectral_response_terms(response::FiniteWidthCollisionLinearization)
    return response.spectral
end
function finite_width_unsupported_offset_terms(response::FiniteWidthCollisionLinearization)
    return response.unsupported_offset
end
function finite_width_unsupported_distribution_terms(
    response::FiniteWidthCollisionLinearization
)
    return response.unsupported_distribution
end

function Base.isempty(response::FiniteWidthCollisionLinearization)
    return isempty(response.terms) &&
           isempty(response.unsupported_offset) &&
           isempty(response.unsupported_distribution)
end

function _append_finite_width_linearization_component!(
    combined::Dict{
        FiniteWidthFrequencyTerm{S,E1,E2},BackgroundOccupationLinearization{D,S}
    },
    occupation_response::Dict{
        FiniteWidthFrequencyTerm{S,E1,E2},BackgroundOccupationLinearization{D,S}
    },
    spectral_response::Dict{
        FiniteWidthFrequencyTerm{S,E1,E2},BackgroundOccupationLinearization{D,S}
    },
    unsupported::Dict{FiniteWidthFrequencyTerm{S,E1,E2},D},
    expression::FiniteWidthFrequencyExpression{C,S,E1,E2,G,Ctx},
    target::FieldFamily{S},
    model::LorentzianSpectralModel{M,S},
    background::OccupationBackground{B},
    response::OccupationSpectralResponse{R,S},
    ω_external::A,
    include_external::Bool,
    ::Type{D},
) where {
    C<:Number,
    S<:Statistics,
    E1,
    E2,
    G,
    Ctx<:AbstractWignerContext,
    M<:Real,
    B<:Number,
    R<:Real,
    A<:Real,
    D<:Number,
}
    for (term, coefficient) in expression
        if finite_width_reduction_kind(term) === FiniteWidthUnsupported
            unsupported[term] = convert(D, coefficient)
            continue
        end

        polynomial = _finite_width_occupation_polynomial(
            term, coefficient, target, include_external
        )
        weight = _finite_width_weight(term, model, ω_external)
        linearization = evaluate_occupation_linearization(
            occupation_linearization(polynomial), background
        )
        occupation_part = _scale_background_linearization(linearization, weight, D)

        background_polynomial = convert(
            D, evaluate_occupation_polynomial(polynomial, background)
        )
        spectral_terms = Pair{OccupationAtom{S},D}[]
        sizehint!(spectral_terms, length(response))
        for (atom, variation) in response
            δweight = _finite_width_weight_variation(term, model, variation, ω_external)
            push!(spectral_terms, atom => background_polynomial * convert(D, δweight))
        end
        spectral_part = BackgroundOccupationLinearization{D,S}(spectral_terms)
        combined_part = _combine_background_linearizations(
            occupation_part, spectral_part, D
        )

        _push_background_linearization!(occupation_response, term, occupation_part, D)
        _push_background_linearization!(spectral_response, term, spectral_part, D)
        _push_background_linearization!(combined, term, combined_part, D)
    end
    return nothing
end

"""
    linearize_finite_width_collision(collision, model, background, spectral_response, ω_external)

Construct the self-consistent finite-width collision Fréchet derivative around `background`.
The fixed-spectral-model occupation derivative and the explicitly supplied spectral-model response
are composed exactly by the product rule.

The spectral response is supplied per occupation variation channel. For a supplied channel, every
spectral line required by a resolved term must have an explicit `(δE, δΓ, δZ)`, including zeros;
missing line data remain an error. The microscopic `ω_external` is held fixed and is not a
collective center-time response frequency. This operation performs no loop-momentum integration,
trap closure, moment projection, or Chapman--Enskog approximation.
"""
function linearize_finite_width_collision(
    collision::FiniteWidthFrequencyCollision{C,S,O,E1,E2,G,Ctx},
    model::LorentzianSpectralModel{M,S},
    background::OccupationBackground{B},
    response::OccupationSpectralResponse{R,S},
    ω_external::A,
) where {
    C<:Number,
    S<:Statistics,
    O,
    E1,
    E2,
    G,
    Ctx<:AbstractWignerContext,
    M<:Real,
    B<:Number,
    R<:Real,
    A<:Real,
}
    D = promote_type(C, M, B, R, A, Rational{Int})
    combined = Dict{
        FiniteWidthFrequencyTerm{S,E1,E2},BackgroundOccupationLinearization{D,S}
    }()
    occupation_response = Dict{
        FiniteWidthFrequencyTerm{S,E1,E2},BackgroundOccupationLinearization{D,S}
    }()
    spectral_response = Dict{
        FiniteWidthFrequencyTerm{S,E1,E2},BackgroundOccupationLinearization{D,S}
    }()
    unsupported_offset = Dict{FiniteWidthFrequencyTerm{S,E1,E2},D}()
    unsupported_distribution = Dict{FiniteWidthFrequencyTerm{S,E1,E2},D}()
    target = target_family(collision)

    _append_finite_width_linearization_component!(
        combined,
        occupation_response,
        spectral_response,
        unsupported_offset,
        collision_offset(collision),
        target,
        model,
        background,
        response,
        ω_external,
        false,
        D,
    )
    _append_finite_width_linearization_component!(
        combined,
        occupation_response,
        spectral_response,
        unsupported_distribution,
        collision_distribution_coefficient(collision),
        target,
        model,
        background,
        response,
        ω_external,
        true,
        D,
    )

    return FiniteWidthCollisionLinearization{D,S,O,E1,E2,G,Ctx}(
        combined,
        occupation_response,
        spectral_response,
        unsupported_offset,
        unsupported_distribution,
        target,
        parameters(collision),
        wigner_context(collision),
    )
end
