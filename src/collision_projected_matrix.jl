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
