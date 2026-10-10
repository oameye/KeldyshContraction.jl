"""First variation of one Lorentzian spectral line's physical data."""
struct LorentzianSpectralDataVariation{C<:Real}
    energy::C
    linewidth::C
    residue::C
end

function LorentzianSpectralDataVariation(
    energy::A, linewidth::B, residue::D
) where {A<:Real,B<:Real,D<:Real}
    C = promote_type(A, B, D)
    return LorentzianSpectralDataVariation{C}(
        convert(C, energy), convert(C, linewidth), convert(C, residue)
    )
end

"""Return the first variation of the Lorentzian line center."""
spectral_energy_variation(variation::LorentzianSpectralDataVariation) = variation.energy

"""Return the first variation of the Lorentzian linewidth."""
function spectral_linewidth_variation(variation::LorentzianSpectralDataVariation)
    return variation.linewidth
end

"""Return the first variation of the Lorentzian residue."""
spectral_residue_variation(variation::LorentzianSpectralDataVariation) = variation.residue

"""Explicit line-by-line first variation of a Lorentzian spectral model."""
struct LorentzianSpectralModelVariation{C<:Real,S<:Statistics}
    data::Dict{SpectralLineIdentity{S},LorentzianSpectralDataVariation{C}}
end

function LorentzianSpectralModelVariation(
    data::AbstractDict{SpectralLineIdentity{S},LorentzianSpectralDataVariation{C}}
) where {C<:Real,S<:Statistics}
    return LorentzianSpectralModelVariation{C,S}(Dict(data))
end

"""Return the explicit spectral-data variation for `line`."""
function spectral_data_variation(
    variation::LorentzianSpectralModelVariation{C,S}, line::SpectralLineIdentity{S}
) where {C<:Real,S<:Statistics}
    haskey(variation.data, line) || throw(
        ArgumentError(
            "finite-width response requires an explicit spectral-data variation for every resolved line",
        ),
    )
    return variation.data[line]
end

"""
Exact first variation of the integrated `m`th power of one Lorentzian line.

For

    B_m = binomial(2m-2,m-1) Z^m / Γ^(m-1),

the derivative is evaluated directly, without dividing by `Z`. The integrated factor is
independent of the line center `E`, so the energy variation does not contribute.
"""
function lorentzian_integrated_power_variation(
    data::LorentzianSpectralData{C},
    variation::LorentzianSpectralDataVariation{V},
    multiplicity::Integer,
) where {C<:Real,V<:Real}
    multiplicity >= 1 || throw(ArgumentError("spectral multiplicity must be positive"))
    m = Int(multiplicity)
    D = promote_type(C, V)
    Z = convert(D, spectral_residue(data))
    Γ = convert(D, spectral_linewidth(data))
    δZ = convert(D, spectral_residue_variation(variation))
    δΓ = convert(D, spectral_linewidth_variation(variation))
    coefficient = convert(D, binomial(2m - 2, m - 1))

    residue_term = coefficient * convert(D, m) * Z^(m - 1) * δZ / Γ^(m - 1)
    m == 1 && return residue_term
    linewidth_term = -coefficient * convert(D, m - 1) * Z^m * δΓ / Γ^m
    return residue_term + linewidth_term
end

function _factorized_spectral_factor_values(
    reduction::LorentzianSpectralReduction{S},
    model::LorentzianSpectralModel{C,S},
    variation::LorentzianSpectralModelVariation{V,S},
) where {C<:Real,V<:Real,S<:Statistics}
    _require_resolved(reduction)
    factors = spectral_factors(reduction)
    D = promote_type(C, V)
    values = Vector{D}(undef, length(factors))
    derivatives = Vector{D}(undef, length(factors))
    for index in eachindex(factors)
        factor = factors[index]
        line = spectral_line(factor)
        data = spectral_data(model, line)
        δdata = spectral_data_variation(variation, line)
        multiplicity = spectral_multiplicity(factor)
        values[index] = convert(D, lorentzian_integrated_power(data, multiplicity))
        derivatives[index] = convert(
            D, lorentzian_integrated_power_variation(data, δdata, multiplicity)
        )
    end
    return values, derivatives
end

"""Exact directional derivative of a resolved factorized Lorentzian spectral weight."""
function evaluate_spectral_weight_variation(
    reduction::LorentzianSpectralReduction{S},
    model::LorentzianSpectralModel{C,S},
    variation::LorentzianSpectralModelVariation{V,S},
) where {C<:Real,V<:Real,S<:Statistics}
    values, derivatives = _factorized_spectral_factor_values(reduction, model, variation)
    D = promote_type(C, V)
    jacobian = convert(D, spectral_jacobian(reduction))
    isempty(values) && return zero(D)

    derivative = zero(D)
    for varied in eachindex(values)
        contribution = derivatives[varied]
        for index in eachindex(values)
            index == varied && continue
            contribution *= values[index]
        end
        derivative += contribution
    end
    return jacobian * derivative
end

function _convolution_residue_and_variation(
    reduction::LorentzianConvolutionReduction{S},
    model::LorentzianSpectralModel{C,S},
    variation::LorentzianSpectralModelVariation{V,S},
) where {C<:Real,V<:Real,S<:Statistics}
    D = promote_type(C, V)
    dependent = convolution_dependent_line(reduction)
    dependent_data = spectral_data(model, dependent)
    dependent_variation = spectral_data_variation(variation, dependent)
    residue = convert(D, spectral_residue(dependent_data))
    δresidue = convert(D, spectral_residue_variation(dependent_variation))

    for line in convolution_pivot_lines(reduction)
        data = spectral_data(model, line)
        δdata = spectral_data_variation(variation, line)
        factor = convert(D, spectral_residue(data))
        δfactor = convert(D, spectral_residue_variation(δdata))
        δresidue = δresidue * factor + residue * δfactor
        residue *= factor
    end
    return residue, δresidue
end

"""
Exact directional derivative of a resolved one-residual-shell Lorentzian convolution.

The external microscopic spectral frequency is held fixed. Variations of line energies,
linewidths and residues are propagated through the exact stored mismatch, effective linewidth and
residue product. This operation does not introduce a center-time response frequency.
"""
function evaluate_spectral_convolution_variation(
    reduction::LorentzianConvolutionReduction{S},
    model::LorentzianSpectralModel{C,S},
    variation::LorentzianSpectralModelVariation{V,S},
    ω_external::A,
) where {A<:Real,C<:Real,V<:Real,S<:Statistics}
    _require_convolution_resolved(reduction)
    D = promote_type(A, C, V)
    dependent = convolution_dependent_line(reduction)
    dependent_data = spectral_data(model, dependent)
    dependent_variation = spectral_data_variation(variation, dependent)

    mismatch =
        convert(D, convolution_external_coefficient(reduction)) * convert(D, ω_external) -
        convert(D, spectral_energy(dependent_data))
    δmismatch = -convert(D, spectral_energy_variation(dependent_variation))
    width = convert(D, spectral_linewidth(dependent_data))
    δwidth = convert(D, spectral_linewidth_variation(dependent_variation))

    coefficients = convolution_coefficients(reduction)
    pivot_lines = convolution_pivot_lines(reduction)
    for index in eachindex(pivot_lines, coefficients)
        line = pivot_lines[index]
        coefficient = coefficients[index]
        data = spectral_data(model, line)
        δdata = spectral_data_variation(variation, line)
        mismatch += convert(D, coefficient) * convert(D, spectral_energy(data))
        δmismatch += convert(D, coefficient) * convert(D, spectral_energy_variation(δdata))
        width += convert(D, abs(coefficient)) * convert(D, spectral_linewidth(data))
        δwidth +=
            convert(D, abs(coefficient)) * convert(D, spectral_linewidth_variation(δdata))
    end

    residue, δresidue = _convolution_residue_and_variation(reduction, model, variation)
    denominator = mismatch^2 + width^2 / convert(D, 4)
    δdenominator = convert(D, 2) * mismatch * δmismatch + width * δwidth / convert(D, 2)
    numerator = residue * width
    δnumerator = δresidue * width + residue * δwidth
    jacobian = convert(D, convolution_jacobian(reduction))
    return jacobian * (δnumerator * denominator - numerator * δdenominator) / denominator^2
end
