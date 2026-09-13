"""Return the number of loop frequencies still active in a partially integrated state."""
function loop_frequency_count(state::FrequencyIntegrationState)
    return length(state.active_loop_basis_indices)
end

function _active_spectral_loop_matrix(state::FrequencyIntegrationState)
    term = state.term
    lines = kinetic_lines(term.carrier)
    kinds = spectral_dispersive_kinds(term)
    spectral_lines = Int[]
    for line_index in state.active_line_indices
        kinds[line_index] === CollisionSpectral && push!(spectral_lines, line_index)
    end

    matrix = zeros(
        EnergyCoefficient, length(spectral_lines), length(state.active_loop_basis_indices)
    )
    for (row, line_index) in enumerate(spectral_lines)
        routed = momentum(lines[line_index])
        for (column, basis_index) in enumerate(state.active_loop_basis_indices)
            matrix[row, column] = routed[basis_index]
        end
    end
    return matrix, spectral_lines
end

"""Rank of the active explicit spectral constraints with respect to active loop frequencies."""
function spectral_frequency_rank(state::FrequencyIntegrationState)::Int
    matrix, _ = _active_spectral_loop_matrix(state)
    return _exact_frequency_rank(matrix)
end

function _require_unshifted_active_lines(state::FrequencyIntegrationState)
    lines = kinetic_lines(state.term.carrier)
    for line_index in state.active_line_indices
        iszero(regularisation_shift(lines[line_index])) || throw(
            ArgumentError(
                "active shifted line requires further Trotter reduction before ordinary frequency reduction",
            ),
        )
    end
    return nothing
end

function _active_full_rank_frequency_reduction(
    state::FrequencyIntegrationState{S}, target::FieldFamily{S}
) where {S<:Statistics}
    _require_unshifted_active_lines(state)
    term = state.term
    lines = kinetic_lines(term.carrier)
    kinds = spectral_dispersive_kinds(term)
    basis = momentum_basis(term)
    external_index = _external_frequency_basis_index(term)
    loop_indices = state.active_loop_basis_indices
    target_momentum = basis_momentum(basis, external_index)
    target_energy = EnergyForm(DispersionAtom(target, target_momentum))

    matrix, spectral_lines = _active_spectral_loop_matrix(state)
    nloops = length(loop_indices)
    rank = _exact_frequency_rank(matrix)
    rank == nloops || throw(
        ArgumentError(
            "active spectral frequency rank $rank does not span $nloops active loop-frequency variables",
        ),
    )

    pivot_rows = _first_full_rank_rows(matrix, nloops)
    pivot_lines = Int[spectral_lines[row] for row in pivot_rows]
    pivot_matrix = matrix[pivot_rows, :]
    inverse, determinant = _exact_inverse_and_determinant(pivot_matrix)

    pivot_rhs = EnergyForm{S}[
        _line_frequency_rhs(lines[line_index], target_energy, external_index) for
        line_index in pivot_lines
    ]
    loop_energies = Vector{EnergyForm{S}}(undef, nloops)
    basis_size = length(basis)
    for loop in 1:nloops
        loop_energies[loop] = _combine_frequency_energies(
            view(inverse, loop, :), pivot_rhs, basis_size
        )
    end

    factor = inv(abs(determinant))
    shells = EnergyShell{S}[]
    principal_values = PrincipalValueSupport{S}[]
    pivot_mask = falses(length(lines))
    for line_index in pivot_lines
        pivot_mask[line_index] = true
    end

    for line_index in state.active_line_indices
        kinds[line_index] === CollisionSpectral && pivot_mask[line_index] && continue
        mismatch = _frequency_mismatch(
            lines[line_index], loop_indices, loop_energies, target_energy, external_index
        )
        if kinds[line_index] === CollisionSpectral
            iszero(mismatch) && throw(
                ArgumentError(
                    "active residual spectral constraint collapsed to δ(0); dependent-shell reduction is required",
                ),
            )
            shell, support_factor = energy_shell(mismatch)
            factor *= support_factor
            push!(shells, shell)
        else
            iszero(mismatch) && throw(
                ArgumentError(
                    "active dispersive denominator collapsed to PV(1/0); a dedicated causal identity is required",
                ),
            )
            support, support_factor = principal_value_support(mismatch)
            factor *= support_factor
            push!(principal_values, support)
        end
    end

    return FullRankFrequencyReduction{S}(
        copy(loop_indices),
        loop_energies,
        pivot_lines,
        FrequencySupport(shells, principal_values),
        factor,
    )
end

function _active_partial_spectral_frequency_reduction(
    state::FrequencyIntegrationState{S}, target::FieldFamily{S}
) where {S<:Statistics}
    _require_unshifted_active_lines(state)
    term = state.term
    lines = kinetic_lines(term.carrier)
    kinds = spectral_dispersive_kinds(term)
    basis = momentum_basis(term)
    external_index = _external_frequency_basis_index(term)
    loop_indices = state.active_loop_basis_indices
    target_momentum = basis_momentum(basis, external_index)
    target_energy = EnergyForm(DispersionAtom(target, target_momentum))

    matrix, spectral_lines = _active_spectral_loop_matrix(state)
    rank = _exact_frequency_rank(matrix)
    pivot_rows = _first_full_rank_rows(matrix, rank)
    pivot_lines = Int[spectral_lines[row] for row in pivot_rows]
    selected_rows = matrix[pivot_rows, :]
    pivot_columns = _first_full_rank_columns(selected_rows, rank)
    residual_columns = _complement_indices(length(loop_indices), pivot_columns)

    pivot_matrix = selected_rows[:, pivot_columns]
    inverse, determinant = _exact_inverse_and_determinant(pivot_matrix)
    residual_matrix = selected_rows[:, residual_columns]
    pivot_residual = _pivot_residual_coefficients(inverse, residual_matrix)

    pivot_rhs = EnergyForm{S}[
        _line_frequency_rhs(lines[line_index], target_energy, external_index) for
        line_index in pivot_lines
    ]
    basis_size = length(basis)
    pivot_offsets = Vector{EnergyForm{S}}(undef, rank)
    for pivot in 1:rank
        pivot_offsets[pivot] = _combine_frequency_energies(
            view(inverse, pivot, :), pivot_rhs, basis_size
        )
    end

    factor = inv(abs(determinant))
    shells = EnergyShell{S}[]
    causal_terms = CausalFrequencyTerm{EnergyCoefficient,S}[CausalFrequencyTerm(
        one(EnergyCoefficient), CausalFrequencyDenominator{S}[]
    )]
    pivot_mask = falses(length(lines))
    for line_index in pivot_lines
        pivot_mask[line_index] = true
    end

    for line_index in state.active_line_indices
        kinds[line_index] === CollisionSpectral && pivot_mask[line_index] && continue
        residual, energy = _partial_frequency_mismatch(
            lines[line_index],
            loop_indices,
            pivot_columns,
            residual_columns,
            pivot_offsets,
            pivot_residual,
            target_energy,
            external_index,
        )

        if kinds[line_index] === CollisionSpectral
            all(iszero, residual) || throw(
                ErrorException(
                    "maximal active spectral elimination unexpectedly left spectral loop-frequency dependence",
                ),
            )
            iszero(energy) && throw(
                ArgumentError(
                    "active residual spectral constraint collapsed to δ(0); dependent-shell reduction is required",
                ),
            )
            shell, support_factor = energy_shell(energy)
            factor *= support_factor
            push!(shells, shell)
        else
            iszero(energy) &&
                all(iszero, residual) &&
                throw(
                    ArgumentError(
                        "active dispersive denominator collapsed to PV(1/0); a dedicated causal identity is required",
                    ),
                )
            causal_terms = _append_dispersive_causal_factor(causal_terms, residual, energy)
        end
    end

    return PartialSpectralFrequencyReduction{S}(
        Int[loop_indices[column] for column in pivot_columns],
        Int[loop_indices[column] for column in residual_columns],
        pivot_lines,
        pivot_offsets,
        pivot_residual,
        FrequencySupport(shells, PrincipalValueSupport{S}[]),
        causal_terms,
        factor,
    )
end

function _active_dependent_spectral_frequency_reduction(
    state::FrequencyIntegrationState{S}, target::FieldFamily{S}
) where {S<:Statistics}
    _require_unshifted_active_lines(state)
    term = state.term
    lines = kinetic_lines(term.carrier)
    kinds = spectral_dispersive_kinds(term)
    basis = momentum_basis(term)
    external_index = _external_frequency_basis_index(term)
    loop_indices = state.active_loop_basis_indices
    target_momentum = basis_momentum(basis, external_index)
    target_energy = EnergyForm(DispersionAtom(target, target_momentum))

    constraints, spectral_lines, _ = _active_spectral_constraints(state, target)
    system = _canonical_spectral_constraint_system(constraints)
    isempty(system.dependent_rows) && throw(
        ArgumentError(
            "active dependent spectral reduction requires linearly dependent shell support",
        ),
    )

    independent_rows = system.independent_rows
    nloops = length(loop_indices)
    loop_matrix = zeros(EnergyCoefficient, length(independent_rows), nloops)
    for (row, system_row) in enumerate(independent_rows)
        constraint = system.constraints[system_row]
        for column in 1:nloops
            loop_matrix[row, column] = constraint.loop_coefficients[column]
        end
    end

    loop_rank = _exact_frequency_rank(loop_matrix)
    pivot_basis_rows = _first_full_rank_rows(loop_matrix, loop_rank)
    pivot_system_rows = Int[independent_rows[row] for row in pivot_basis_rows]
    selected_rows = loop_matrix[pivot_basis_rows, :]
    pivot_columns = _first_full_rank_columns(selected_rows, loop_rank)
    residual_columns = _complement_indices(nloops, pivot_columns)

    pivot_matrix = selected_rows[:, pivot_columns]
    inverse, determinant = _exact_inverse_and_determinant(pivot_matrix)
    residual_matrix = selected_rows[:, residual_columns]
    pivot_residual = _pivot_residual_coefficients(inverse, residual_matrix)

    pivot_rhs = EnergyForm{S}[system.constraints[row].rhs for row in pivot_system_rows]
    basis_size = length(basis)
    pivot_offsets = Vector{EnergyForm{S}}(undef, loop_rank)
    for pivot in 1:loop_rank
        pivot_offsets[pivot] = _combine_frequency_energies(
            view(inverse, pivot, :), pivot_rhs, basis_size
        )
    end

    factor = _constraint_system_scale_factor(system) * inv(abs(determinant))
    shells = EnergyShell{S}[]
    causal_terms = CausalFrequencyTerm{EnergyCoefficient,S}[CausalFrequencyTerm(
        one(EnergyCoefficient), CausalFrequencyDenominator{S}[]
    )]

    pivot_system_mask = falses(length(system.constraints))
    for row in pivot_system_rows
        pivot_system_mask[row] = true
    end
    for system_row in independent_rows
        pivot_system_mask[system_row] && continue
        residual, energy = _canonical_constraint_frequency_mismatch(
            system.constraints[system_row],
            pivot_columns,
            residual_columns,
            pivot_offsets,
            pivot_residual,
        )
        all(iszero, residual) || throw(
            ErrorException(
                "affine-independent active spectral quotient unexpectedly retained loop-frequency dependence",
            ),
        )
        iszero(energy) && throw(
            ErrorException(
                "affine-independent active spectral row collapsed after frequency elimination",
            ),
        )
        shell, support_factor = energy_shell(energy)
        factor *= support_factor
        push!(shells, shell)
    end

    for line_index in state.active_line_indices
        kinds[line_index] === CollisionDispersive || continue
        residual, energy = _partial_frequency_mismatch(
            lines[line_index],
            loop_indices,
            pivot_columns,
            residual_columns,
            pivot_offsets,
            pivot_residual,
            target_energy,
            external_index,
        )
        iszero(energy) &&
            all(iszero, residual) &&
            throw(
                ArgumentError(
                    "active dispersive denominator collapsed to PV(1/0); a dedicated causal identity is required",
                ),
            )
        causal_terms = _append_dispersive_causal_factor(causal_terms, residual, energy)
    end

    pivot_lines = Int[spectral_lines[system.source_rows[row]] for row in pivot_system_rows]
    return PartialSpectralFrequencyReduction{S}(
        Int[loop_indices[column] for column in pivot_columns],
        Int[loop_indices[column] for column in residual_columns],
        pivot_lines,
        pivot_offsets,
        pivot_residual,
        FrequencySupport(shells, PrincipalValueSupport{S}[]),
        causal_terms,
        factor,
    )
end

function _scale_frequency_support(
    reduced::Dict{FrequencySupport{S},ComplexRationals}, factor::ComplexRationals
) where {S<:Statistics}
    isone(factor) && return reduced
    out = Dict{FrequencySupport{S},ComplexRationals}()
    for (support, coefficient) in reduced
        _push_frequency_support_coefficient!(out, support, factor * coefficient)
    end
    return out
end

"""Reduce all remaining unshifted active frequencies after earlier exact integrations."""
function general_frequency_reduction(
    state::FrequencyIntegrationState{S}, target::FieldFamily{S}
) where {S<:Statistics}
    reduction = if spectral_frequency_rank(state) == loop_frequency_count(state)
        full = _active_full_rank_frequency_reduction(state, target)
        Dict(frequency_support(full) => convert(ComplexRationals, frequency_factor(full)))
    else
        reduce_partial_spectral_frequency(
            _active_partial_spectral_frequency_reduction(state, target)
        )
    end
    return _scale_frequency_support(reduction, state.factor)
end

"""Reduce the finite quotient multiplying active dependent-shell support."""
function dependent_spectral_frequency_reduction(
    state::FrequencyIntegrationState{S}, target::FieldFamily{S}
) where {S<:Statistics}
    reduction = reduce_partial_spectral_frequency(
        _active_dependent_spectral_frequency_reduction(state, target)
    )
    return _scale_frequency_support(reduction, state.factor)
end
