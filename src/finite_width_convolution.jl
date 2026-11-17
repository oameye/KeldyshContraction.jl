"""Structural outcome of a one-residual-shell Lorentzian convolution reduction."""
@enum LorentzianConvolutionReductionKind::UInt8 begin
    LorentzianConvolutionResolved = 0
    LorentzianConvolutionUnsupported = 1
end

"""
Concrete analytic reduction of `L + 1` distinct spectral lines over `L` loop frequencies.

For a resolved reduction, `pivot_lines` define the independent Lorentzian frequency coordinates,
`dependent_lines` contains the single residual line, and `coefficients` gives the exact affine
relation between the dependent frequency and the pivot frequencies. `external_coefficient`
multiplies the projected external frequency. `mismatch` is exactly the residual quasiparticle
energy mismatch used by the strict frequency reducer.
"""
struct LorentzianConvolutionReduction{S<:Statistics}
    kind::LorentzianConvolutionReductionKind
    jacobian::EnergyCoefficient
    pivot_lines::Vector{SpectralLineIdentity{S}}
    dependent_lines::Vector{SpectralLineIdentity{S}}
    coefficients::Vector{EnergyCoefficient}
    external_coefficient::EnergyCoefficient
    mismatch::EnergyForm{S}
end

convolution_reduction_kind(reduction::LorentzianConvolutionReduction) = reduction.kind
function convolution_reduction_resolved(reduction::LorentzianConvolutionReduction)
    return reduction.kind === LorentzianConvolutionResolved
end

function _unsupported_lorentzian_convolution(
    ::Type{S}, basis_size::Integer
) where {S<:Statistics}
    return LorentzianConvolutionReduction{S}(
        LorentzianConvolutionUnsupported,
        zero(EnergyCoefficient),
        SpectralLineIdentity{S}[],
        SpectralLineIdentity{S}[],
        EnergyCoefficient[],
        zero(EnergyCoefficient),
        EnergyForm{S}(basis_size),
    )
end

function _require_convolution_resolved(reduction::LorentzianConvolutionReduction)
    convolution_reduction_resolved(reduction) || throw(
        ArgumentError(
            "Lorentzian convolution reduction is unsupported for this frequency structure",
        ),
    )
    return nothing
end

function convolution_jacobian(reduction::LorentzianConvolutionReduction)
    _require_convolution_resolved(reduction)
    return reduction.jacobian
end
function convolution_pivot_lines(reduction::LorentzianConvolutionReduction)
    _require_convolution_resolved(reduction)
    return reduction.pivot_lines
end
function convolution_dependent_line(reduction::LorentzianConvolutionReduction)
    _require_convolution_resolved(reduction)
    return only(reduction.dependent_lines)
end
function convolution_coefficients(reduction::LorentzianConvolutionReduction)
    _require_convolution_resolved(reduction)
    return reduction.coefficients
end
function convolution_external_coefficient(reduction::LorentzianConvolutionReduction)
    _require_convolution_resolved(reduction)
    return reduction.external_coefficient
end
function convolution_energy_mismatch(reduction::LorentzianConvolutionReduction)
    _require_convolution_resolved(reduction)
    return reduction.mismatch
end

"""
    lorentzian_convolution_reduction(term, target)

Resolve the analytic Lorentzian convolution when an all-spectral term has exactly one more
distinct spectral-line identity than independent loop frequencies. The line frequencies are
reduced with the same exact routing matrix used by the strict quasiparticle reducer.

The operation is structural and deterministic. Repeated lines, dispersive factors, finite
Trotter shifts, higher-codimension residual shells, and rank-deficient systems return a concrete
`LorentzianConvolutionUnsupported` result. No numerical integration or regulator is introduced.
"""
function lorentzian_convolution_reduction(
    term::SpectralDispersiveTerm{S}, target::FieldFamily{S}
) where {S<:Statistics}
    basis = momentum_basis(term)
    unsupported = _unsupported_lorentzian_convolution(S, length(basis))
    lines = kinetic_lines(term.carrier)
    kinds = spectral_dispersive_kinds(term)
    all(kind -> kind === CollisionSpectral, kinds) || return unsupported
    any(line -> !iszero(regularisation_shift(line)), lines) && return unsupported

    identities = SpectralLineIdentity{S}[_spectral_line_identity(line) for line in lines]
    length(unique(identities)) == length(identities) || return unsupported

    matrix, spectral_lines = _spectral_loop_matrix(term)
    loop_indices = _loop_frequency_basis_indices(term)
    nloops = length(loop_indices)
    length(spectral_lines) == nloops + 1 || return unsupported
    _exact_frequency_rank(matrix) == nloops || return unsupported

    pivot_rows = _first_full_rank_rows(matrix, nloops)
    pivot_mask = falses(length(spectral_lines))
    for row in pivot_rows
        pivot_mask[row] = true
    end
    dependent_rows = Int[row for row in eachindex(spectral_lines) if !pivot_mask[row]]
    length(dependent_rows) == 1 || return unsupported
    dependent_row = only(dependent_rows)

    pivot_matrix = matrix[pivot_rows, :]
    inverse, determinant = _exact_inverse_and_determinant(pivot_matrix)
    coefficients = Vector{EnergyCoefficient}(undef, nloops)
    for pivot in 1:nloops
        coefficient = zero(EnergyCoefficient)
        for loop in 1:nloops
            coefficient += matrix[dependent_row, loop] * inverse[loop, pivot]
        end
        coefficients[pivot] = coefficient
    end

    external_index = _external_frequency_basis_index(term)
    pivot_lines = Int[spectral_lines[row] for row in pivot_rows]
    dependent_line = spectral_lines[dependent_row]
    external_coefficient = momentum(lines[dependent_line])[external_index]
    for pivot in eachindex(pivot_lines)
        external_coefficient -=
            coefficients[pivot] * momentum(lines[pivot_lines[pivot]])[external_index]
    end

    target_momentum = basis_momentum(basis, external_index)
    target_energy = EnergyForm(DispersionAtom(target, target_momentum))
    pivot_rhs = EnergyForm{S}[
        _line_frequency_rhs(lines[line_index], target_energy, external_index) for
        line_index in pivot_lines
    ]
    loop_energies = Vector{EnergyForm{S}}(undef, nloops)
    for loop in 1:nloops
        loop_energies[loop] = _combine_frequency_energies(
            view(inverse, loop, :), pivot_rhs, length(basis)
        )
    end
    mismatch = _frequency_mismatch(
        lines[dependent_line], loop_indices, loop_energies, target_energy, external_index
    )
    iszero(mismatch) && return unsupported

    return LorentzianConvolutionReduction{S}(
        LorentzianConvolutionResolved,
        inv(abs(determinant)),
        SpectralLineIdentity{S}[_spectral_line_identity(lines[i]) for i in pivot_lines],
        SpectralLineIdentity{S}[_spectral_line_identity(lines[dependent_line])],
        coefficients,
        external_coefficient,
        mismatch,
    )
end

"""
Return the exact effective linewidth of a resolved Lorentzian convolution.

If the dependent frequency is `ω_d = βω_ext + Σ c_j ω_j`, Cauchy stability gives
`Γ_eff = Γ_d + Σ |c_j| Γ_j`.
"""
function convolution_effective_linewidth(
    reduction::LorentzianConvolutionReduction{S}, model::LorentzianSpectralModel{C,S}
) where {C<:Real,S<:Statistics}
    _require_convolution_resolved(reduction)
    width = spectral_linewidth(model, only(reduction.dependent_lines))
    for i in eachindex(reduction.pivot_lines, reduction.coefficients)
        coefficient = abs(reduction.coefficients[i])
        width +=
            convert(C, coefficient) * spectral_linewidth(model, reduction.pivot_lines[i])
    end
    return width
end

"""
Evaluate a resolved one-residual-shell Lorentzian convolution at external energy `ω_external`.

The exact result is the affine frequency Jacobian times one Lorentzian in the residual energy
mismatch, with linewidth `Γ_eff = Γ_d + Σ |c_j|Γ_j` and residue equal to the product of the line
residues. Missing spectral data or an unsupported reduction is an error.
"""
function evaluate_spectral_convolution(
    reduction::LorentzianConvolutionReduction{S},
    model::LorentzianSpectralModel{C,S},
    ω_external::A,
) where {A<:Real,C<:Real,S<:Statistics}
    _require_convolution_resolved(reduction)
    D = promote_type(A, C)
    dependent = only(reduction.dependent_lines)
    dependent_data = spectral_data(model, dependent)

    mismatch =
        convert(D, reduction.external_coefficient) * convert(D, ω_external) -
        convert(D, spectral_energy(dependent_data))
    width = convert(D, spectral_linewidth(dependent_data))
    residue = convert(D, spectral_residue(dependent_data))
    for i in eachindex(reduction.pivot_lines, reduction.coefficients)
        data = spectral_data(model, reduction.pivot_lines[i])
        coefficient = reduction.coefficients[i]
        mismatch += convert(D, coefficient) * convert(D, spectral_energy(data))
        width += convert(D, abs(coefficient)) * convert(D, spectral_linewidth(data))
        residue *= convert(D, spectral_residue(data))
    end

    return convert(D, reduction.jacobian) * residue * width /
           (mismatch^2 + width^2 / convert(D, 4))
end
