function _finite_width_external_spectral_line(
    term::FiniteWidthFrequencyTerm{S}, target::FieldFamily{S}
) where {S<:Statistics}
    source = finite_width_source_term(term)
    basis = momentum_basis(source)
    external_index = _external_frequency_basis_index(source)
    external_momentum = basis_momentum(basis, external_index)
    return SpectralLineIdentity{S}(target, external_momentum)
end

function _external_projected_finite_width_weight(
    term::FiniteWidthFrequencyTerm{S},
    model::LorentzianSpectralModel{M,S},
    target::FieldFamily{S},
) where {S<:Statistics,M<:Real}
    kind = finite_width_reduction_kind(term)
    external_line = _finite_width_external_spectral_line(term, target)
    if kind === FiniteWidthFactorized
        internal = evaluate_spectral_weight(finite_width_spectral_reduction(term), model)
        return spectral_residue(model, external_line) * internal
    elseif kind === FiniteWidthConvolution
        return evaluate_external_spectral_projection(
            finite_width_convolution_reduction(term), model, external_line
        )
    end
    return throw(
        ArgumentError("unsupported finite-width term has no external spectral weight")
    )
end

function _external_projected_finite_width_weight_variation(
    term::FiniteWidthFrequencyTerm{S},
    model::LorentzianSpectralModel{M,S},
    variation::LorentzianSpectralModelVariation{V,S},
    target::FieldFamily{S},
) where {S<:Statistics,M<:Real,V<:Real}
    kind = finite_width_reduction_kind(term)
    external_line = _finite_width_external_spectral_line(term, target)
    if kind === FiniteWidthFactorized
        internal = evaluate_spectral_weight(finite_width_spectral_reduction(term), model)
        δinternal = evaluate_spectral_weight_variation(
            finite_width_spectral_reduction(term), model, variation
        )
        external_data = spectral_data(model, external_line)
        external_variation = spectral_data_variation(variation, external_line)
        Z = spectral_residue(external_data)
        δZ = spectral_residue_variation(external_variation)
        return δZ * internal + Z * δinternal
    elseif kind === FiniteWidthConvolution
        return evaluate_external_spectral_projection_variation(
            finite_width_convolution_reduction(term), model, variation, external_line
        )
    end
    return throw(
        ArgumentError(
            "unsupported finite-width term has no external spectral-weight variation"
        ),
    )
end

"""
Self-consistent collision response after integrating the microscopic external spectral line.

For every analytically resolved finite-width term, `terms` stores the complete product-rule
response, `occupation` stores the fixed-spectral occupation contribution, and `spectral` stores
the response of the full spectral measure including the external line. Unsupported frequency
structures remain explicit.

This representation is intentionally distinct from [`FiniteWidthCollisionLinearization`](@ref),
which is evaluated at a fixed microscopic external frequency.
"""
struct ExternalSpectralCollisionLinearization{
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

order(::ExternalSpectralCollisionLinearization{C,S,O}) where {C,S,O} = O
statistics(::ExternalSpectralCollisionLinearization{C,S}) where {C,S} = S
parameters(response::ExternalSpectralCollisionLinearization) = response.parameter
target_family(response::ExternalSpectralCollisionLinearization) = response.target
function gradient_order(
    ::ExternalSpectralCollisionLinearization{C,S,O,E1,E2,G}
) where {C,S,O,E1,E2,G}
    return Val(G)
end
wigner_context(response::ExternalSpectralCollisionLinearization) = response.context

"""Return the complete response after microscopic external spectral projection."""
function external_spectral_linearized_terms(
    response::ExternalSpectralCollisionLinearization
)
    return response.terms
end

"""Return the fixed-spectral occupation part after external spectral projection."""
function external_spectral_occupation_response_terms(
    response::ExternalSpectralCollisionLinearization
)
    return response.occupation
end

"""Return the full spectral-model response, including the external line."""
function external_spectral_model_response_terms(
    response::ExternalSpectralCollisionLinearization
)
    return response.spectral
end

function finite_width_unsupported_offset_terms(
    response::ExternalSpectralCollisionLinearization
)
    return response.unsupported_offset
end
function finite_width_unsupported_distribution_terms(
    response::ExternalSpectralCollisionLinearization
)
    return response.unsupported_distribution
end

function Base.isempty(response::ExternalSpectralCollisionLinearization)
    return isempty(response.terms) &&
           isempty(response.unsupported_offset) &&
           isempty(response.unsupported_distribution)
end

function _append_external_spectral_linearization_component!(
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
        weight = _external_projected_finite_width_weight(term, model, target)
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
            δweight = _external_projected_finite_width_weight_variation(
                term, model, variation, target
            )
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
    linearize_external_spectral_collision(collision, model, background, spectral_response)

Construct the exact occupation-space Fréchet derivative after analytically integrating the
remaining microscopic external spectral line.

For a convolution sector this uses the full Cauchy width
`Γ_int + |β|Γ_ext`; for a factorized sector the normalized external spectral integral supplies
its residue `Z_ext`. The spectral response differentiates the external line together with every
internal line. The result therefore represents the full microscopic spectral measure rather than
a collision evaluated at fixed `ω_external`.

This operation performs no loop-momentum integration, trap closure, moment projection,
Chapman--Enskog approximation, or collective center-time-frequency construction.
"""
function linearize_external_spectral_collision(
    collision::FiniteWidthFrequencyCollision{C,S,O,E1,E2,G,Ctx},
    model::LorentzianSpectralModel{M,S},
    background::OccupationBackground{B},
    response::OccupationSpectralResponse{R,S},
) where {
    C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext,M<:Real,B<:Number,R<:Real
}
    D = promote_type(C, M, B, R, Rational{Int})
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

    _append_external_spectral_linearization_component!(
        combined,
        occupation_response,
        spectral_response,
        unsupported_offset,
        collision_offset(collision),
        target,
        model,
        background,
        response,
        false,
        D,
    )
    _append_external_spectral_linearization_component!(
        combined,
        occupation_response,
        spectral_response,
        unsupported_distribution,
        collision_distribution_coefficient(collision),
        target,
        model,
        background,
        response,
        true,
        D,
    )

    return ExternalSpectralCollisionLinearization{D,S,O,E1,E2,G,Ctx}(
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
