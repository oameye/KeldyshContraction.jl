"""
Exact partial elimination of independent spectral frequency constraints.

The selected spectral delta functions eliminate a maximal set of loop frequencies. Remaining
loop-frequency dependence is carried only by exact causal denominators. This representation is
independent of perturbative order, graph topology, and line count.
"""
struct PartialSpectralFrequencyReduction{S<:Statistics}
    pivot_loop_basis_indices::Vector{Int}
    residual_loop_basis_indices::Vector{Int}
    pivot_spectral_lines::Vector{Int}
    pivot_frequency_offsets::Vector{EnergyForm{S}}
    pivot_residual_coefficients::Matrix{EnergyCoefficient}
    support::FrequencySupport{S}
    causal_terms::Vector{CausalFrequencyTerm{EnergyCoefficient,S}}
    factor::EnergyCoefficient
end

function Base.isequal(
    a::PartialSpectralFrequencyReduction{S}, b::PartialSpectralFrequencyReduction{S}
) where {S<:Statistics}
    return isequal(a.pivot_loop_basis_indices, b.pivot_loop_basis_indices) &&
           isequal(a.residual_loop_basis_indices, b.residual_loop_basis_indices) &&
           isequal(a.pivot_spectral_lines, b.pivot_spectral_lines) &&
           isequal(a.pivot_frequency_offsets, b.pivot_frequency_offsets) &&
           isequal(a.pivot_residual_coefficients, b.pivot_residual_coefficients) &&
           isequal(a.support, b.support) &&
           isequal(a.causal_terms, b.causal_terms) &&
           a.factor == b.factor
end
function Base.:(==)(
    a::PartialSpectralFrequencyReduction, b::PartialSpectralFrequencyReduction
)
    return isequal(a, b)
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
    return throw(
        ArgumentError("spectral constraints do not contain the requested column rank")
    )
end

function _complement_indices(length_total::Int, selected::Vector{Int})
    mask = falses(length_total)
    for index in selected
        mask[index] = true
    end
    return Int[index for index in 1:length_total if !mask[index]]
end

function _pivot_residual_coefficients(
    inverse::Matrix{EnergyCoefficient}, residual_matrix::Matrix{EnergyCoefficient}
)
    rank = size(inverse, 1)
    nresidual = size(residual_matrix, 2)
    out = zeros(EnergyCoefficient, rank, nresidual)
    for row in 1:rank
        for column in 1:nresidual
            value = zero(EnergyCoefficient)
            for inner in 1:rank
                value -= inverse[row, inner] * residual_matrix[inner, column]
            end
            out[row, column] = value
        end
    end
    return out
end

function _partial_frequency_mismatch(
    line::KineticLine{S},
    loop_indices::Vector{Int},
    pivot_columns::Vector{Int},
    residual_columns::Vector{Int},
    pivot_offsets::Vector{EnergyForm{S}},
    pivot_residual::Matrix{EnergyCoefficient},
    target_energy::EnergyForm{S},
    external_index::Int,
) where {S<:Statistics}
    routed = momentum(line)
    energy =
        routed[external_index] * target_energy -
        EnergyForm(DispersionAtom(line.family, routed))
    residual = EnergyCoefficient[
        routed[loop_indices[column]] for column in residual_columns
    ]

    for pivot in eachindex(pivot_columns)
        coefficient = routed[loop_indices[pivot_columns[pivot]]]
        iszero(coefficient) && continue
        energy = energy + coefficient * pivot_offsets[pivot]
        for residual_index in eachindex(residual)
            residual[residual_index] += coefficient * pivot_residual[pivot, residual_index]
        end
    end
    return residual, energy
end

function _canonical_constraint_frequency_mismatch(
    constraint::AffineSpectralConstraint{S},
    pivot_columns::Vector{Int},
    residual_columns::Vector{Int},
    pivot_offsets::Vector{EnergyForm{S}},
    pivot_residual::Matrix{EnergyCoefficient},
) where {S<:Statistics}
    energy = -constraint.rhs
    residual = EnergyCoefficient[
        constraint.loop_coefficients[column] for column in residual_columns
    ]
    for pivot in eachindex(pivot_columns)
        coefficient = constraint.loop_coefficients[pivot_columns[pivot]]
        iszero(coefficient) && continue
        energy = energy + coefficient * pivot_offsets[pivot]
        for residual_index in eachindex(residual)
            residual[residual_index] += coefficient * pivot_residual[pivot, residual_index]
        end
    end
    return residual, energy
end

function _append_dispersive_causal_factor(
    terms::Vector{CausalFrequencyTerm{EnergyCoefficient,S}},
    loop_coefficients::Vector{EnergyCoefficient},
    energy::EnergyForm{S},
) where {S<:Statistics}
    next = CausalFrequencyTerm{EnergyCoefficient,S}[]
    sizehint!(next, 2 * length(terms))
    for term in terms
        for infinitesimal in (-1 // 1, 1 // 1)
            denominators = copy(term.denominators)
            push!(
                denominators,
                CausalFrequencyDenominator{S}(
                    copy(loop_coefficients), energy, infinitesimal
                ),
            )
            push!(next, CausalFrequencyTerm(term.coefficient / 2, denominators))
        end
    end
    return next
end

"""
    partial_spectral_frequency_reduction(term, target)

Use a deterministic maximal independent subset of the explicit spectral constraints to
eliminate as many loop frequencies as possible. The remaining dispersive factors are expanded
as exact retarded/advanced causal denominators while retaining their infinitesimal
prescriptions. No residual loop-frequency integral is performed by this function.

Nonzero equal-time shifts remain outside this generic path and require the dedicated Trotter
rule.
"""
function partial_spectral_frequency_reduction(
    term::SpectralDispersiveTerm{S}, target::FieldFamily{S}
) where {S<:Statistics}
    lines = kinetic_lines(term.carrier)
    for line in lines
        iszero(regularisation_shift(line)) || throw(
            ArgumentError(
                "shifted equal-time line requires dedicated Trotter frequency reduction"
            ),
        )
    end

    basis = momentum_basis(term)
    external_index = _external_frequency_basis_index(term)
    loop_indices = _loop_frequency_basis_indices(term)
    target_momentum = basis_momentum(basis, external_index)
    target_energy = EnergyForm(DispersionAtom(target, target_momentum))

    matrix, spectral_lines = _spectral_loop_matrix(term)
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

    kinds = spectral_dispersive_kinds(term)
    for line_index in eachindex(lines)
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
                    "maximal spectral elimination unexpectedly left spectral loop-frequency dependence",
                ),
            )
            iszero(energy) && throw(
                ArgumentError(
                    "residual spectral constraint collapsed to δ(0); a dedicated repeated-shell rule is required",
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
                        "dispersive denominator collapsed to PV(1/0); a dedicated causal identity is required",
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

"""
    dependent_spectral_frequency_reduction(term, target)

Reduce an unshifted term with linearly dependent spectral support in two exact stages. First,
canonicalize the complete affine shell system and remove dependent rows in the quotient of
`(loop frequencies, symbolic energies)`. Frequency pivots are then chosen only from that
independent affine basis. This prevents a redundant shell from reappearing as a duplicated
residual `EnergyShell` under a different pivot choice.

The returned finite factor contains the product of every original delta-row scale, the exact
frequency-pivot Jacobian, residual independent shell normalisations, and subsequent causal
factors. The singular relation itself remains in `DependentShellSupport` and is never assigned
a value.
"""
function dependent_spectral_frequency_reduction(
    term::SpectralDispersiveTerm{S}, target::FieldFamily{S}
) where {S<:Statistics}
    lines = kinetic_lines(term.carrier)
    for line in lines
        iszero(regularisation_shift(line)) || throw(
            ArgumentError(
                "shifted equal-time line requires dedicated Trotter frequency reduction"
            ),
        )
    end

    basis = momentum_basis(term)
    external_index = _external_frequency_basis_index(term)
    loop_indices = _loop_frequency_basis_indices(term)
    target_momentum = basis_momentum(basis, external_index)
    target_energy = EnergyForm(DispersionAtom(target, target_momentum))

    constraints, spectral_lines, _ = _term_spectral_constraints(term, target)
    system = _canonical_spectral_constraint_system(constraints)
    isempty(system.dependent_rows) && throw(
        ArgumentError(
            "dependent spectral reduction requires linearly dependent shell support"
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
                "affine-independent spectral quotient unexpectedly retained loop-frequency dependence",
            ),
        )
        iszero(energy) && throw(
            ErrorException(
                "affine-independent spectral row collapsed after frequency elimination"
            ),
        )
        shell, support_factor = energy_shell(energy)
        factor *= support_factor
        push!(shells, shell)
    end

    kinds = spectral_dispersive_kinds(term)
    for line_index in eachindex(lines)
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
                    "dispersive denominator collapsed to PV(1/0); a dedicated causal identity is required",
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

function _combine_frequency_support(
    a::FrequencySupport{S}, b::FrequencySupport{S}
) where {S}
    return FrequencySupport(
        vcat(a.shells, b.shells), vcat(a.principal_values, b.principal_values)
    )
end

function _complex_causal_term(
    term::CausalFrequencyTerm{C,S}
) where {C<:Number,S<:Statistics}
    return CausalFrequencyTerm(
        convert(ComplexRationals, term.coefficient), term.denominators
    )
end

"""Integrate all residual causal loop frequencies and lower the result to canonical support."""
function reduce_partial_spectral_frequency(
    reduction::PartialSpectralFrequencyReduction{S}
) where {S<:Statistics}
    terms = CausalFrequencyTerm{ComplexRationals,S}[
        _complex_causal_term(term) for term in reduction.causal_terms
    ]
    nresidual = length(reduction.residual_loop_basis_indices)
    for frequency_index in 1:nresidual
        next = CausalFrequencyTerm{ComplexRationals,S}[]
        for term in terms
            append!(next, integrate_causal_frequency(term, frequency_index))
        end
        terms = next
        isempty(terms) && break
    end

    out = Dict{FrequencySupport{S},ComplexRationals}()
    prefactor = convert(ComplexRationals, reduction.factor)
    for term in terms
        lowered = lower_causal_frequency_term(term)
        for (support, coefficient) in lowered
            combined_support = _combine_frequency_support(reduction.support, support)
            _push_frequency_support_coefficient!(
                out, combined_support, prefactor * coefficient
            )
        end
    end
    return out
end

"""
Reduce one unshifted spectral/dispersive term through the exact available frequency machinery.

Full explicit spectral rank uses the existing #290 fast path unchanged. Rank-deficient terms
first eliminate their maximal spectral subset and then use causal residue integration for the
remaining loop-frequency subspace.
"""
function general_frequency_reduction(
    term::SpectralDispersiveTerm{S}, target::FieldFamily{S}
) where {S<:Statistics}
    if spectral_frequency_rank(term) == loop_frequency_count(term)
        reduction = full_rank_frequency_reduction(term, target)
        return Dict(
            frequency_support(reduction) =>
                convert(ComplexRationals, frequency_factor(reduction)),
        )
    end
    return reduce_partial_spectral_frequency(
        partial_spectral_frequency_reduction(term, target)
    )
end
