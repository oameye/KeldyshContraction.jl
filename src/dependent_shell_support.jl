"""
One exact affine spectral constraint

    c⋅ω = E

before loop-frequency elimination. The loop coefficients and symbolic energy right-hand side are
kept separately so the same object can serve both integration and singular-support analysis.
"""
struct AffineSpectralConstraint{S<:Statistics}
    loop_coefficients::Vector{EnergyCoefficient}
    rhs::EnergyForm{S}
end

function AffineSpectralConstraint(
    loop_coefficients::AbstractVector{<:Rational}, rhs::EnergyForm{S}
) where {S<:Statistics}
    return AffineSpectralConstraint{S}(
        EnergyCoefficient[convert(EnergyCoefficient, value) for value in loop_coefficients],
        rhs,
    )
end

function Base.isequal(
    a::AffineSpectralConstraint{S}, b::AffineSpectralConstraint{S}
) where {S<:Statistics}
    return isequal(a.loop_coefficients, b.loop_coefficients) && isequal(a.rhs, b.rhs)
end
Base.:(==)(a::AffineSpectralConstraint, b::AffineSpectralConstraint) = isequal(a, b)
function Base.hash(constraint::AffineSpectralConstraint, h::UInt)
    return hash(
        constraint.rhs,
        hash(constraint.loop_coefficients, hash(AffineSpectralConstraint, h)),
    )
end

function _constraint_nonzero_count(constraint::AffineSpectralConstraint)
    return count(value -> !iszero(value), constraint.loop_coefficients) +
           length(constraint.rhs)
end

function Base.isless(
    a::AffineSpectralConstraint{S}, b::AffineSpectralConstraint{S}
) where {S<:Statistics}
    a_nnz = _constraint_nonzero_count(a)
    b_nnz = _constraint_nonzero_count(b)
    a_nnz == b_nnz || return a_nnz < b_nnz
    length(a.loop_coefficients) == length(b.loop_coefficients) ||
        return length(a.loop_coefficients) < length(b.loop_coefficients)
    for index in eachindex(a.loop_coefficients)
        a_value = a.loop_coefficients[index]
        b_value = b.loop_coefficients[index]
        a_value == b_value || return a_value < b_value
    end
    return isless(a.rhs, b.rhs)
end

function _constraint_scale(constraint::AffineSpectralConstraint)
    for coefficient in constraint.loop_coefficients
        iszero(coefficient) || return coefficient
    end
    iszero(constraint.rhs) && return one(EnergyCoefficient)
    return last(first(energy_terms(constraint.rhs)))
end

function _canonical_spectral_constraint(constraint::AffineSpectralConstraint{S}) where {S}
    scale = _constraint_scale(constraint)
    inverse_scale = inv(scale)
    canonical = AffineSpectralConstraint{S}(
        EnergyCoefficient[inverse_scale * value for value in constraint.loop_coefficients],
        inverse_scale * constraint.rhs,
    )
    return canonical, inv(abs(scale))
end

"""
Canonical physical support of a linearly dependent spectral-constraint system.

`independent_constraints` is a deterministic exact basis of the affine shell system. Each row
of `dependency_rows` gives the rational coefficients expressing one additional canonical
constraint as a linear combination of that basis.

Scale/Jacobian factors are deliberately not stored here: rescaling a delta constraint changes
the coefficient multiplying the support, not the support identity itself.
"""
struct DependentShellSupport{S<:Statistics}
    independent_constraints::Vector{AffineSpectralConstraint{S}}
    dependency_rows::Vector{Vector{EnergyCoefficient}}
end

function Base.isequal(a::DependentShellSupport{S}, b::DependentShellSupport{S}) where {S}
    return isequal(a.independent_constraints, b.independent_constraints) &&
           isequal(a.dependency_rows, b.dependency_rows)
end
Base.:(==)(a::DependentShellSupport, b::DependentShellSupport) = isequal(a, b)
function Base.hash(support::DependentShellSupport, h::UInt)
    return hash(
        support.dependency_rows,
        hash(support.independent_constraints, hash(DependentShellSupport, h)),
    )
end

constraint_rank(support::DependentShellSupport) = length(support.independent_constraints)
function constraint_count(support::DependentShellSupport)
    return constraint_rank(support) + length(support.dependency_rows)
end
function has_dependent_shell_support(support::DependentShellSupport)
    return !isempty(support.dependency_rows)
end

function repeated_shell_multiplicity(support::DependentShellSupport)
    constraint_rank(support) == 1 || return 0
    return constraint_count(support)
end

"""
Canonical dependent-shell support together with the exact factor generated solely by
canonicalising the individual delta constraints.

`constraint_scale` does not include any loop-frequency integration Jacobian or residual
shell/PV normalisation. Those belong to the subsequent frequency-reduction stage.
"""
struct DependentShellReduction{S<:Statistics}
    support::DependentShellSupport{S}
    constraint_scale::EnergyCoefficient
end

dependent_shell_support(reduction::DependentShellReduction) = reduction.support
constraint_scale_factor(reduction::DependentShellReduction) = reduction.constraint_scale
function has_dependent_shell_support(reduction::DependentShellReduction)
    return has_dependent_shell_support(reduction.support)
end
constraint_rank(reduction::DependentShellReduction) = constraint_rank(reduction.support)
constraint_count(reduction::DependentShellReduction) = constraint_count(reduction.support)
function repeated_shell_multiplicity(reduction::DependentShellReduction)
    return repeated_shell_multiplicity(reduction.support)
end

function _spectral_constraint_atoms(
    constraints::Vector{AffineSpectralConstraint{S}}
) where {S<:Statistics}
    atoms = DispersionAtom{S}[]
    for constraint in constraints
        for (atom, _) in constraint.rhs
            push!(atoms, atom)
        end
    end
    sort!(atoms)
    unique!(atoms)
    return atoms
end

function _spectral_constraint_matrix(
    constraints::Vector{AffineSpectralConstraint{S}}
) where {S<:Statistics}
    isempty(constraints) && return zeros(EnergyCoefficient, 0, 0)
    nloops = length(first(constraints).loop_coefficients)
    basis_size = energy_basis_size(first(constraints).rhs)
    for constraint in constraints
        length(constraint.loop_coefficients) == nloops ||
            throw(DimensionMismatch("spectral constraints use different loop dimensions"))
        energy_basis_size(constraint.rhs) == basis_size ||
            throw(DimensionMismatch("spectral constraints use different energy bases"))
    end

    atoms = _spectral_constraint_atoms(constraints)
    atom_index = Dict{DispersionAtom{S},Int}(
        atom => index for (index, atom) in pairs(atoms)
    )
    matrix = zeros(EnergyCoefficient, length(constraints), nloops + length(atoms))
    for row in eachindex(constraints)
        constraint = constraints[row]
        for column in 1:nloops
            matrix[row, column] = constraint.loop_coefficients[column]
        end
        for (atom, coefficient) in constraint.rhs
            matrix[row, nloops + atom_index[atom]] = -coefficient
        end
    end
    return matrix
end

"""
Canonicalized affine spectral system with source-row provenance.

Rows are sorted by physical canonical constraint before exact rank analysis. `source_rows`
therefore maps each canonical row back to the corresponding row of the caller's input. The
independent/dependent partition is computed in the *full affine space* containing both loop-
frequency coefficients and symbolic energy atoms; later frequency pivots must be chosen only
inside `independent_rows`.
"""
struct CanonicalSpectralConstraintSystem{S<:Statistics}
    constraints::Vector{AffineSpectralConstraint{S}}
    source_rows::Vector{Int}
    scale_factors::Vector{EnergyCoefficient}
    matrix::Matrix{EnergyCoefficient}
    independent_rows::Vector{Int}
    dependent_rows::Vector{Int}
end

function _constraint_order_isless(
    i::Int, j::Int, constraints::Vector{AffineSpectralConstraint{S}}
) where {S<:Statistics}
    a = constraints[i]
    b = constraints[j]
    isequal(a, b) && return i < j
    return isless(a, b)
end

function _canonical_spectral_constraint_system(
    input::AbstractVector{AffineSpectralConstraint{S}}
) where {S<:Statistics}
    canonical_unsorted = AffineSpectralConstraint{S}[]
    scale_unsorted = EnergyCoefficient[]
    sizehint!(canonical_unsorted, length(input))
    sizehint!(scale_unsorted, length(input))
    for constraint in input
        canonical, scale_factor = _canonical_spectral_constraint(constraint)
        push!(canonical_unsorted, canonical)
        push!(scale_unsorted, scale_factor)
    end

    order = collect(eachindex(canonical_unsorted))
    sort!(order; lt=(i, j) -> _constraint_order_isless(i, j, canonical_unsorted))
    constraints = AffineSpectralConstraint{S}[canonical_unsorted[row] for row in order]
    scale_factors = EnergyCoefficient[scale_unsorted[row] for row in order]
    source_rows = Int[row for row in order]
    matrix = _spectral_constraint_matrix(constraints)
    rank = _exact_frequency_rank(matrix)
    independent_rows = _first_full_rank_rows(matrix, rank)
    independent_mask = falses(length(constraints))
    for row in independent_rows
        independent_mask[row] = true
    end
    dependent_rows = Int[row for row in eachindex(constraints) if !independent_mask[row]]
    return CanonicalSpectralConstraintSystem{S}(
        constraints, source_rows, scale_factors, matrix, independent_rows, dependent_rows
    )
end

function _constraint_system_scale_factor(system::CanonicalSpectralConstraintSystem)
    factor = one(EnergyCoefficient)
    for scale_factor in system.scale_factors
        factor *= scale_factor
    end
    return factor
end

function _constraint_dependency_row(
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
            throw(ErrorException("exact spectral dependency reconstruction failed"))
    end
    return coefficients
end

function _dependent_shell_reduction(system::CanonicalSpectralConstraintSystem{S}) where {S}
    independent = AffineSpectralConstraint{S}[
        system.constraints[row] for row in system.independent_rows
    ]
    dependencies = Vector{EnergyCoefficient}[]
    rank = length(system.independent_rows)
    if rank == 0
        for _ in system.dependent_rows
            push!(dependencies, EnergyCoefficient[])
        end
        support = DependentShellSupport{S}(independent, dependencies)
        return DependentShellReduction{S}(support, _constraint_system_scale_factor(system))
    end

    pivot_matrix = system.matrix[system.independent_rows, :]
    pivot_columns = _first_full_rank_columns(pivot_matrix, rank)
    square = pivot_matrix[:, pivot_columns]
    inverse, _ = _exact_inverse_and_determinant(square)
    for row in system.dependent_rows
        push!(
            dependencies,
            _constraint_dependency_row(
                view(system.matrix, row, :), pivot_matrix, pivot_columns, inverse
            ),
        )
    end
    support = DependentShellSupport{S}(independent, dependencies)
    return DependentShellReduction{S}(support, _constraint_system_scale_factor(system))
end

"""
    analyze_spectral_dependencies(constraints)

Canonicalise an arbitrary exact spectral-constraint system and preserve every rational linear
dependency instead of turning it into a `δ(0)` value. The returned support is invariant under
input row ordering and under nonzero rational rescaling/sign reversal of individual
constraints. `constraint_scale_factor` returns the product of the individual delta scaling
factors only; it deliberately excludes frequency-integration and residual-support Jacobians.
"""
function analyze_spectral_dependencies(
    input::AbstractVector{AffineSpectralConstraint{S}}
) where {S<:Statistics}
    return _dependent_shell_reduction(_canonical_spectral_constraint_system(input))
end

"""Concrete dependency analysis of the spectral lines in one collision term."""
struct SpectralDependencyAnalysis{S<:Statistics}
    reduction::DependentShellReduction{S}
    spectral_line_indices::Vector{Int}
    regularisation_shifts::Vector{Int8}
end

function dependent_shell_support(analysis::SpectralDependencyAnalysis)
    return dependent_shell_support(analysis.reduction)
end
function constraint_scale_factor(analysis::SpectralDependencyAnalysis)
    return constraint_scale_factor(analysis.reduction)
end
function has_dependent_shell_support(analysis::SpectralDependencyAnalysis)
    return has_dependent_shell_support(analysis.reduction)
end
constraint_rank(analysis::SpectralDependencyAnalysis) = constraint_rank(analysis.reduction)
function constraint_count(analysis::SpectralDependencyAnalysis)
    return constraint_count(analysis.reduction)
end
function repeated_shell_multiplicity(analysis::SpectralDependencyAnalysis)
    return repeated_shell_multiplicity(analysis.reduction)
end

function _term_spectral_constraints(
    term::SpectralDispersiveTerm{S}, target::FieldFamily{S}
) where {S<:Statistics}
    lines = kinetic_lines(term.carrier)
    basis = momentum_basis(term)
    external_index = _external_frequency_basis_index(term)
    loop_indices = _loop_frequency_basis_indices(term)
    target_momentum = basis_momentum(basis, external_index)
    target_energy = EnergyForm(DispersionAtom(target, target_momentum))
    kinds = spectral_dispersive_kinds(term)

    constraints = AffineSpectralConstraint{S}[]
    line_indices = Int[]
    shifts = Int8[]
    for line_index in eachindex(lines)
        kinds[line_index] === CollisionSpectral || continue
        line = lines[line_index]
        coefficients = EnergyCoefficient[momentum(line)[index] for index in loop_indices]
        rhs = _line_frequency_rhs(line, target_energy, external_index)
        push!(constraints, AffineSpectralConstraint{S}(coefficients, rhs))
        push!(line_indices, line_index)
        push!(shifts, regularisation_shift(line))
    end
    return constraints, line_indices, shifts
end

"""
Analyse exact linear dependence of all explicit spectral constraints in a collision term.

Equal-time/Trotter shifts are retained as separate provenance and do not participate in the
mass-shell dependency relation.
"""
function analyze_spectral_dependencies(
    term::SpectralDispersiveTerm{S}, target::FieldFamily{S}
) where {S<:Statistics}
    constraints, line_indices, shifts = _term_spectral_constraints(term, target)
    return SpectralDependencyAnalysis{S}(
        analyze_spectral_dependencies(constraints), line_indices, shifts
    )
end
