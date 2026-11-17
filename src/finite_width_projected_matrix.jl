function _require_finite_width_projection_resolved(response)
    isempty(finite_width_unsupported_offset_terms(response)) || throw(
        ArgumentError(
            "cannot project finite-width response with unresolved offset sectors"
        ),
    )
    isempty(finite_width_unsupported_distribution_terms(response)) || throw(
        ArgumentError(
            "cannot project finite-width response with unresolved distribution sectors"
        ),
    )
    return nothing
end

function _finite_width_projection_terms(response::FiniteWidthCollisionLinearization)
    return finite_width_linearized_terms(response)
end
function _finite_width_projection_terms(response::ExternalSpectralCollisionLinearization)
    return external_spectral_linearized_terms(response)
end

function _finite_width_projection_resolved(response)
    _require_finite_width_projection_resolved(response)
    return _finite_width_projection_terms(response)
end

function _projected_finite_width_entry(
    terms,
    left,
    right,
    closure::CollisionPerturbationClosure,
    functional::CollisionProjectionFunctional{V},
)::V where {V<:Number}
    entry = zero(V)
    for (sector, linearization) in terms
        for (atom, coefficient) in linearization
            amplitude = perturbation_amplitude(closure, right, atom)
            iszero(amplitude) && continue
            closed_response = coefficient * amplitude
            entry += project_collision_channel(
                functional, left, sector, atom, closed_response
            )
        end
    end
    return entry
end

function _projected_finite_width_row(
    terms,
    left,
    right::Tuple,
    closure::CollisionPerturbationClosure,
    functional::CollisionProjectionFunctional,
)
    return map(
        descriptor ->
            _projected_finite_width_entry(terms, left, descriptor, closure, functional),
        right,
    )
end

function _projected_finite_width_matrix(
    response,
    basis::CollisionProjectionBasis{L,R},
    closure::CollisionPerturbationClosure,
    functional::CollisionProjectionFunctional{V},
) where {L<:Tuple,R<:Tuple,V<:Number}
    terms = _finite_width_projection_resolved(response)
    left = left_projection_basis(basis)
    right = right_perturbation_basis(basis)
    rows = map(
        descriptor ->
            _projected_finite_width_row(terms, descriptor, right, closure, functional),
        left,
    )
    values = Matrix{V}(undef, length(left), length(right))
    for i in eachindex(rows), j in eachindex(right)
        values[i, j] = rows[i][j]
    end
    return ProjectedCollisionMatrix{V,L,R}(left, right, values)
end

"""
    projected_collision_matrix(response, basis, closure, functional)

Assemble the finite-dimensional projection of a fully resolved self-consistent finite-width
collision response.

Both fixed-microscopic-frequency [`FiniteWidthCollisionLinearization`](@ref) and fully externally
projected [`ExternalSpectralCollisionLinearization`](@ref) use the same response-agnostic basis,
right perturbation closure, and left projection functional. The supplied response already contains
its exact occupation and spectral-data chain rule.

No trap model, quadrature, Chapman--Enskog closure, collective frequency, or pole condition is
hidden here. Projection is rejected if unresolved finite-width offset or distribution sectors
remain, so unsupported collision physics cannot be silently discarded.
"""
function projected_collision_matrix(
    response::FiniteWidthCollisionLinearization,
    basis::CollisionProjectionBasis,
    closure::CollisionPerturbationClosure,
    functional::CollisionProjectionFunctional,
)
    return _projected_finite_width_matrix(response, basis, closure, functional)
end

function projected_collision_matrix(
    response::ExternalSpectralCollisionLinearization,
    basis::CollisionProjectionBasis,
    closure::CollisionPerturbationClosure,
    functional::CollisionProjectionFunctional,
)
    return _projected_finite_width_matrix(response, basis, closure, functional)
end
