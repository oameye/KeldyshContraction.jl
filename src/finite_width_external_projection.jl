function _external_spectral_projection_geometry(
    reduction::LorentzianConvolutionReduction{S},
    model::LorentzianSpectralModel{C,S},
    external_line::SpectralLineIdentity{S},
) where {C<:Real,S<:Statistics}
    _require_convolution_resolved(reduction)
    external_data = spectral_data(model, external_line)
    dependent = convolution_dependent_line(reduction)
    dependent_data = spectral_data(model, dependent)
    β = convolution_external_coefficient(reduction)

    mismatch =
        convert(C, β) * spectral_energy(external_data) - spectral_energy(dependent_data)
    width =
        convert(C, abs(β)) * spectral_linewidth(external_data) +
        spectral_linewidth(dependent_data)
    residue = spectral_residue(external_data) * spectral_residue(dependent_data)

    coefficients = convolution_coefficients(reduction)
    pivot_lines = convolution_pivot_lines(reduction)
    for index in eachindex(pivot_lines, coefficients)
        data = spectral_data(model, pivot_lines[index])
        coefficient = coefficients[index]
        mismatch += convert(C, coefficient) * spectral_energy(data)
        width += convert(C, abs(coefficient)) * spectral_linewidth(data)
        residue *= spectral_residue(data)
    end
    return mismatch, width, residue
end

"""
    external_spectral_projection_mismatch(reduction, model, external_line)

Return the residual energy mismatch after analytically projecting the remaining microscopic
external frequency against `external_line`.
"""
function external_spectral_projection_mismatch(
    reduction::LorentzianConvolutionReduction,
    model::LorentzianSpectralModel,
    external_line::SpectralLineIdentity,
)
    mismatch, _, _ = _external_spectral_projection_geometry(reduction, model, external_line)
    return mismatch
end

"""
    external_spectral_projection_linewidth(reduction, model, external_line)

Return the full Cauchy linewidth after the external spectral projection. If the internal residual
frequency is `βω_ext + Δ₀`, the result is

    Γ_full = Γ_int + |β| Γ_ext.
"""
function external_spectral_projection_linewidth(
    reduction::LorentzianConvolutionReduction,
    model::LorentzianSpectralModel,
    external_line::SpectralLineIdentity,
)
    _, width, _ = _external_spectral_projection_geometry(reduction, model, external_line)
    return width
end

"""
    evaluate_external_spectral_projection(reduction, model, external_line)

Analytically integrate a resolved internal Lorentzian convolution against the remaining external
spectral line,

    ∫ dω/(2π) A_ext(ω) W_int(ω).

Cauchy stability gives one Lorentzian with the full linewidth
`Γ_int + |β|Γ_ext`, the full energy mismatch, and the product of all line residues. This is a
microscopic spectral projection; it does not introduce a collective center-time response
frequency or any trap/moment closure.
"""
function evaluate_external_spectral_projection(
    reduction::LorentzianConvolutionReduction{S},
    model::LorentzianSpectralModel{C,S},
    external_line::SpectralLineIdentity{S},
) where {C<:Real,S<:Statistics}
    mismatch, width, residue = _external_spectral_projection_geometry(
        reduction, model, external_line
    )
    denominator = mismatch^2 + width^2 / convert(C, 4)
    return convert(C, convolution_jacobian(reduction)) * residue * width / denominator
end

function _external_spectral_projection_geometry_variation(
    reduction::LorentzianConvolutionReduction{S},
    model::LorentzianSpectralModel{C,S},
    variation::LorentzianSpectralModelVariation{V,S},
    external_line::SpectralLineIdentity{S},
) where {C<:Real,V<:Real,S<:Statistics}
    _require_convolution_resolved(reduction)
    D = promote_type(C, V)
    external_data = spectral_data(model, external_line)
    external_variation = spectral_data_variation(variation, external_line)
    dependent = convolution_dependent_line(reduction)
    dependent_data = spectral_data(model, dependent)
    dependent_variation = spectral_data_variation(variation, dependent)
    β = convolution_external_coefficient(reduction)
    βD = convert(D, β)
    absβD = convert(D, abs(β))

    mismatch =
        βD * convert(D, spectral_energy(external_data)) -
        convert(D, spectral_energy(dependent_data))
    δmismatch =
        βD * convert(D, spectral_energy_variation(external_variation)) -
        convert(D, spectral_energy_variation(dependent_variation))
    width =
        absβD * convert(D, spectral_linewidth(external_data)) +
        convert(D, spectral_linewidth(dependent_data))
    δwidth =
        absβD * convert(D, spectral_linewidth_variation(external_variation)) +
        convert(D, spectral_linewidth_variation(dependent_variation))

    external_residue = convert(D, spectral_residue(external_data))
    dependent_residue = convert(D, spectral_residue(dependent_data))
    δexternal_residue = convert(D, spectral_residue_variation(external_variation))
    δdependent_residue = convert(D, spectral_residue_variation(dependent_variation))
    residue = external_residue * dependent_residue
    δresidue = δexternal_residue * dependent_residue + external_residue * δdependent_residue

    coefficients = convolution_coefficients(reduction)
    pivot_lines = convolution_pivot_lines(reduction)
    for index in eachindex(pivot_lines, coefficients)
        line = pivot_lines[index]
        coefficient = coefficients[index]
        data = spectral_data(model, line)
        δdata = spectral_data_variation(variation, line)
        coefficientD = convert(D, coefficient)
        abscoefficientD = convert(D, abs(coefficient))
        mismatch += coefficientD * convert(D, spectral_energy(data))
        δmismatch += coefficientD * convert(D, spectral_energy_variation(δdata))
        width += abscoefficientD * convert(D, spectral_linewidth(data))
        δwidth += abscoefficientD * convert(D, spectral_linewidth_variation(δdata))

        factor = convert(D, spectral_residue(data))
        δfactor = convert(D, spectral_residue_variation(δdata))
        δresidue = δresidue * factor + residue * δfactor
        residue *= factor
    end
    return mismatch, δmismatch, width, δwidth, residue, δresidue
end

"""
    evaluate_external_spectral_projection_variation(reduction, model, variation, external_line)

Exact directional derivative of the full external spectral projection. Variations of line
centers, linewidths, and residues are propagated through the full four-line mismatch, linewidth,
and residue product. The microscopic external frequency has already been integrated out.
"""
function evaluate_external_spectral_projection_variation(
    reduction::LorentzianConvolutionReduction{S},
    model::LorentzianSpectralModel{C,S},
    variation::LorentzianSpectralModelVariation{V,S},
    external_line::SpectralLineIdentity{S},
) where {C<:Real,V<:Real,S<:Statistics}
    mismatch, δmismatch, width, δwidth, residue, δresidue = _external_spectral_projection_geometry_variation(
        reduction, model, variation, external_line
    )
    D = promote_type(C, V)
    denominator = mismatch^2 + width^2 / convert(D, 4)
    δdenominator = convert(D, 2) * mismatch * δmismatch + width * δwidth / convert(D, 2)
    numerator = residue * width
    δnumerator = δresidue * width + residue * δwidth
    jacobian = convert(D, convolution_jacobian(reduction))
    return jacobian * (δnumerator * denominator - numerator * δdenominator) / denominator^2
end
