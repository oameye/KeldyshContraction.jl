"""
One real affine frequency constraint

    a⋅ω + E = 0

underlying a causal denominator. The infinitesimal prescription is deliberately absent: it
selects the boundary value, while the physical singular support is the real affine geometry.
"""
struct AffineFrequencyConstraint{S<:Statistics}
    loop_coefficients::Vector{EnergyCoefficient}
    energy::EnergyForm{S}
end

function AffineFrequencyConstraint(
    loop_coefficients::AbstractVector{<:Rational}, energy::EnergyForm{S}
) where {S<:Statistics}
    return AffineFrequencyConstraint{S}(
        EnergyCoefficient[convert(EnergyCoefficient, value) for value in loop_coefficients],
        energy,
    )
end

function AffineFrequencyConstraint(denominator::CausalFrequencyDenominator{S}) where {S}
    return AffineFrequencyConstraint{S}(
        copy(denominator.loop_coefficients), denominator.energy
    )
end

function Base.isequal(
    a::AffineFrequencyConstraint{S}, b::AffineFrequencyConstraint{S}
) where {S<:Statistics}
    return isequal(a.loop_coefficients, b.loop_coefficients) && isequal(a.energy, b.energy)
end
Base.:(==)(a::AffineFrequencyConstraint, b::AffineFrequencyConstraint) = isequal(a, b)
function Base.hash(constraint::AffineFrequencyConstraint, h::UInt)
    return hash(
        constraint.energy,
        hash(constraint.loop_coefficients, hash(AffineFrequencyConstraint, h)),
    )
end

function Base.isless(
    a::AffineFrequencyConstraint{S}, b::AffineFrequencyConstraint{S}
) where {S<:Statistics}
    length(a.loop_coefficients) == length(b.loop_coefficients) ||
        return length(a.loop_coefficients) < length(b.loop_coefficients)
    for index in eachindex(a.loop_coefficients, b.loop_coefficients)
        a_value = a.loop_coefficients[index]
        b_value = b.loop_coefficients[index]
        a_value == b_value || return a_value < b_value
    end
    return isless(a.energy, b.energy)
end

function _affine_constraint_scale(constraint::AffineFrequencyConstraint)
    for coefficient in constraint.loop_coefficients
        iszero(coefficient) || return coefficient
    end
    for (_, coefficient) in constraint.energy
        iszero(coefficient) || return coefficient
    end
    return one(EnergyCoefficient)
end

function _canonical_affine_constraint(
    constraint::AffineFrequencyConstraint{S}
) where {S<:Statistics}
    scale = _affine_constraint_scale(constraint)
    inverse_scale = inv(scale)
    return AffineFrequencyConstraint{S}(
        EnergyCoefficient[inverse_scale * value for value in constraint.loop_coefficients],
        inverse_scale * constraint.energy,
    )
end

function _affine_constraint_atoms(
    constraints::Vector{AffineFrequencyConstraint{S}}
) where {S<:Statistics}
    atoms = DispersionAtom{S}[]
    for constraint in constraints
        for (atom, _) in constraint.energy
            push!(atoms, atom)
        end
    end
    sort!(atoms)
    unique!(atoms)
    return atoms
end

function _affine_constraint_matrix(
    constraints::Vector{AffineFrequencyConstraint{S}}
) where {S<:Statistics}
    isempty(constraints) && return zeros(EnergyCoefficient, 0, 0)
    nfrequencies = length(first(constraints).loop_coefficients)
    basis_size = energy_basis_size(first(constraints).energy)
    for constraint in constraints
        length(constraint.loop_coefficients) == nfrequencies || throw(
            DimensionMismatch("affine constraints use different frequency dimensions")
        )
        energy_basis_size(constraint.energy) == basis_size ||
            throw(DimensionMismatch("affine constraints use different energy bases"))
    end

    atoms = _affine_constraint_atoms(constraints)
    atom_index = Dict{DispersionAtom{S},Int}(
        atom => index for (index, atom) in pairs(atoms)
    )
    matrix = zeros(EnergyCoefficient, length(constraints), nfrequencies + length(atoms))
    for row in eachindex(constraints)
        constraint = constraints[row]
        for column in 1:nfrequencies
            matrix[row, column] = constraint.loop_coefficients[column]
        end
        for (atom, coefficient) in constraint.energy
            matrix[row, nfrequencies + atom_index[atom]] = coefficient
        end
    end
    return matrix
end

function _first_full_rank_columns(matrix::Matrix{EnergyCoefficient}, target_rank::Int)
    target_rank == 0 && return Int[]
    selected = Int[]
    current_rank = 0
    for column in axes(matrix, 2)
        candidate_columns = vcat(selected, column)
        candidate = matrix[:, candidate_columns]
        candidate_rank = _exact_frequency_rank(candidate)
        if candidate_rank > current_rank
            push!(selected, column)
            current_rank = candidate_rank
            current_rank == target_rank && return selected
        end
    end
    return throw(ArgumentError("affine constraint basis does not have requested rank"))
end

function _affine_dependency_row(
    row::AbstractVector{EnergyCoefficient},
    pivot_matrix::Matrix{EnergyCoefficient},
    pivot_columns::Vector{Int},
    inverse::Matrix{EnergyCoefficient},
)
    rank = size(pivot_matrix, 1)
    coefficients = zeros(EnergyCoefficient, rank)
    for output in 1:rank
        for inner in 1:rank
            coefficients[output] += row[pivot_columns[inner]] * inverse[inner, output]
        end
    end
    for column in axes(pivot_matrix, 2)
        reconstructed = zero(EnergyCoefficient)
        for pivot in 1:rank
            reconstructed += coefficients[pivot] * pivot_matrix[pivot, column]
        end
        reconstructed == row[column] ||
            throw(ErrorException("exact affine dependency reconstruction failed"))
    end
    return coefficients
end

"""
Canonical dependent-active real affine subsystem.

Only rows that participate in at least one exact dependency are retained. Equivalently, a row
belongs to the support iff deleting it leaves the exact row rank unchanged. Independent
spectator denominators are therefore excluded automatically. `dependency_rows` expresses each
additional canonical active constraint in the deterministic basis `independent_constraints`.
"""
struct AffineSingularSupport{S<:Statistics}
    independent_constraints::Vector{AffineFrequencyConstraint{S}}
    dependency_rows::Vector{Vector{EnergyCoefficient}}
end

function Base.isequal(a::AffineSingularSupport{S}, b::AffineSingularSupport{S}) where {S}
    return isequal(a.independent_constraints, b.independent_constraints) &&
           isequal(a.dependency_rows, b.dependency_rows)
end
Base.:(==)(a::AffineSingularSupport, b::AffineSingularSupport) = isequal(a, b)
function Base.hash(support::AffineSingularSupport, h::UInt)
    return hash(
        support.dependency_rows,
        hash(support.independent_constraints, hash(AffineSingularSupport, h)),
    )
end

function affine_support_rank(support::AffineSingularSupport)
    return length(support.independent_constraints)
end
function affine_support_count(support::AffineSingularSupport)
    return affine_support_rank(support) + length(support.dependency_rows)
end
function has_affine_singular_support(support::AffineSingularSupport)
    return !isempty(support.dependency_rows)
end

function _dependent_active_rows(matrix::Matrix{EnergyCoefficient})
    nrows = size(matrix, 1)
    nrows == 0 && return Int[]
    rank = _exact_frequency_rank(matrix)
    active = Int[]
    sizehint!(active, nrows)
    for row in axes(matrix, 1)
        remaining = Int[index for index in axes(matrix, 1) if index != row]
        reduced = matrix[remaining, :]
        _exact_frequency_rank(reduced) == rank && push!(active, row)
    end
    return active
end

function _canonical_affine_constraint_data(
    input::AbstractVector{AffineFrequencyConstraint{S}}
) where {S<:Statistics}
    canonical = AffineFrequencyConstraint{S}[
        _canonical_affine_constraint(constraint) for constraint in input
    ]
    sort!(canonical)
    matrix = _affine_constraint_matrix(canonical)
    return canonical, matrix, _dependent_active_rows(matrix)
end

"""
    affine_singular_support(constraints)

Return the pivot-independent exact affine geometry carried by the dependent-active rows of a
constraint system. Input ordering and nonzero rational rescaling/sign reversal of individual
constraints do not change the returned support.
"""
function affine_singular_support(
    input::AbstractVector{AffineFrequencyConstraint{S}}
) where {S<:Statistics}
    canonical, matrix, active_rows = _canonical_affine_constraint_data(input)
    isempty(active_rows) && return AffineSingularSupport{S}(
        AffineFrequencyConstraint{S}[], Vector{EnergyCoefficient}[]
    )

    active_constraints = AffineFrequencyConstraint{S}[canonical[row] for row in active_rows]
    active_matrix = matrix[active_rows, :]
    rank = _exact_frequency_rank(active_matrix)
    independent_rows = _first_full_rank_rows(active_matrix, rank)
    independent_mask = falses(length(active_rows))
    for row in independent_rows
        independent_mask[row] = true
    end
    dependent_rows = Int[row for row in eachindex(active_rows) if !independent_mask[row]]
    independent = AffineFrequencyConstraint{S}[
        active_constraints[row] for row in independent_rows
    ]

    dependencies = Vector{EnergyCoefficient}[]
    if rank == 0
        for _ in dependent_rows
            push!(dependencies, EnergyCoefficient[])
        end
        return AffineSingularSupport{S}(independent, dependencies)
    end

    pivot_matrix = active_matrix[independent_rows, :]
    pivot_columns = _first_full_rank_columns(pivot_matrix, rank)
    square = pivot_matrix[:, pivot_columns]
    inverse, _ = _exact_inverse_and_determinant(square)
    for row in dependent_rows
        push!(
            dependencies,
            _affine_dependency_row(
                view(active_matrix, row, :), pivot_matrix, pivot_columns, inverse
            ),
        )
    end
    return AffineSingularSupport{S}(independent, dependencies)
end

function affine_singular_support(
    denominators::AbstractVector{CausalFrequencyDenominator{S}}
) where {S<:Statistics}
    constraints = AffineFrequencyConstraint{S}[
        AffineFrequencyConstraint(denominator) for denominator in denominators
    ]
    return affine_singular_support(constraints)
end

function affine_singular_support(term::CausalFrequencyTerm{C,S}) where {C,S<:Statistics}
    return affine_singular_support(causal_frequency_denominators(term))
end
