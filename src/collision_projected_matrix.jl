"""
Explicit finite-dimensional projection basis for a background-linearized collision operator.

The left tuple stores observable/test-function descriptors and the right tuple stores perturbation
basis descriptors. Tuples keep heterogeneous user descriptors concrete and inspectable without
assuming that the two bases coincide.
"""
struct CollisionProjectionBasis{L<:Tuple,R<:Tuple}
    left::L
    right::R
end

"""Return the explicit left test-function basis."""
left_projection_basis(basis::CollisionProjectionBasis) = basis.left

"""Return the explicit right perturbation basis."""
right_perturbation_basis(basis::CollisionProjectionBasis) = basis.right

"""
Typed map from one right-basis descriptor and occupation variation channel to its amplitude.

The evaluator is deliberately user supplied: the collision compiler does not assume a trap model,
phase-space closure, equilibrium ansatz, or numerical quadrature.
"""
struct CollisionPerturbationClosure{V<:Number,F}
    evaluator::F
end

function CollisionPerturbationClosure(evaluator::F, ::Type{V}) where {V<:Number,F}
    return CollisionPerturbationClosure{V,F}(evaluator)
end

"""Evaluate one right-basis perturbation on an exact occupation variation channel."""
function perturbation_amplitude(
    closure::CollisionPerturbationClosure{V}, right, atom::OccupationAtom
)::V where {V<:Number}
    return convert(V, closure.evaluator(right, atom))
end

"""
Typed left projection functional for one already-closed collision-response channel.

Its evaluator receives `(left, sector, atom, response)`, where `response` is the exact
background response coefficient multiplied by the right-basis amplitude. The functional may then
perform an analytic phase-space integral, moment closure, or another explicit downstream
projection. No such rule is hidden in this type.
"""
struct CollisionProjectionFunctional{V<:Number,F}
    evaluator::F
end

function CollisionProjectionFunctional(evaluator::F, ::Type{V}) where {V<:Number,F}
    return CollisionProjectionFunctional{V,F}(evaluator)
end

"""Apply the explicit left projection functional to one closed response channel."""
function project_collision_channel(
    functional::CollisionProjectionFunctional{V},
    left,
    sector::CollisionKernelSector,
    atom::OccupationAtom,
    response,
)::V where {V<:Number}
    return convert(V, functional.evaluator(left, sector, atom, response))
end

"""
Finite-dimensional projected collision matrix with explicit left and right bases.

Rows correspond to `left_projection_basis`, columns to `right_perturbation_basis`. The matrix may
contain numerical or symbolic scalar values. This representation makes no assumption that a later
finite-width projected operator is frequency independent.
"""
struct ProjectedCollisionMatrix{V<:Number,L<:Tuple,R<:Tuple}
    left::L
    right::R
    values::Matrix{V}
end

left_projection_basis(projected::ProjectedCollisionMatrix) = projected.left
right_perturbation_basis(projected::ProjectedCollisionMatrix) = projected.right
matrix(projected::ProjectedCollisionMatrix) = projected.values

function _projected_collision_entry(
    kernel::BackgroundLinearizedCollisionKernel,
    left,
    right,
    closure::CollisionPerturbationClosure,
    functional::CollisionProjectionFunctional{V},
)::V where {V<:Number}
    entry = zero(V)
    for (sector, linearization) in background_linearized_collision_terms(kernel)
        for (atom, coefficient) in linearization
            amplitude = perturbation_amplitude(closure, right, atom)
            iszero(amplitude) && continue
            response = coefficient * amplitude
            entry += project_collision_channel(functional, left, sector, atom, response)
        end
    end
    return entry
end

function _projected_collision_row(
    kernel::BackgroundLinearizedCollisionKernel,
    left,
    right::Tuple,
    closure::CollisionPerturbationClosure,
    functional::CollisionProjectionFunctional,
)
    return map(
        descriptor ->
            _projected_collision_entry(kernel, left, descriptor, closure, functional),
        right,
    )
end

"""
    projected_collision_matrix(kernel, basis, closure, functional)

Assemble the finite-dimensional projection of a supplied-background collision response.

For a background-linearized kernel `δC = Σ_{s,a} c_{s,a} δn_a`, a right-basis closure provides
`ψ_j(a)` and the left functional evaluates the projected response. The resulting entries are
`K_ij = Σ_{s,a} P_i[s,a,c_{s,a} ψ_j(a)]`.

The phase-space projection is therefore explicit and user controlled rather than encoded in the
collision IR.
"""
function projected_collision_matrix(
    kernel::BackgroundLinearizedCollisionKernel{C,S},
    basis::CollisionProjectionBasis{L,R},
    closure::CollisionPerturbationClosure{A},
    functional::CollisionProjectionFunctional{V},
) where {C<:Number,S<:Statistics,L<:Tuple,R<:Tuple,A<:Number,V<:Number}
    left = left_projection_basis(basis)
    right = right_perturbation_basis(basis)
    rows = map(
        descriptor ->
            _projected_collision_row(kernel, descriptor, right, closure, functional),
        left,
    )
    values = Matrix{V}(undef, length(left), length(right))
    for i in eachindex(rows), j in eachindex(right)
        values[i, j] = rows[i][j]
    end
    return ProjectedCollisionMatrix{V,L,R}(left, right, values)
end
