@doc """
    field_families(x)

Return the physical field families carried by an interaction.
""" field_families

@doc """
    target_family(x)

Return the physical field family selected for the external two-point function or kinetic
result.
""" target_family

@doc """
    parameters(x)

Return the perturbative parameter monomial, or the available monomials for a multi-process
result.
""" parameters

@doc """
    matrix(x)

Assemble the statistics-specific Keldysh matrix from the semantic retarded, advanced, and
Keldysh components of a propagator or self-energy.
""" matrix

@doc """
    order(x)

Return the perturbative order represented by `x`.
""" order

@doc """
    statistics(x)

Return the field statistics represented by `x`.
""" statistics

@doc """
    gradient_order(x)

Return the retained Wigner-gradient order as a `Val`. The current kinetic compiler supports
only the homogeneous zeroth-gradient collision problem.
""" gradient_order

@doc """
    wigner_context(x)

Return the Wigner context carried by a Wigner-space or kinetic result.
""" wigner_context

@doc """
    topologies(diagrams)

Group diagrams by their canonical uncolored topology signature. This is an advanced inspection
helper; physical diagram identity still retains field, Keldysh, regularisation, and propagator
information beyond the topology label.
""" topologies

@doc """
    wigner_transform(x; gradient_order=Val(0))

Transform a Fourier-space propagator or self-energy to homogeneous Wigner variables while
preserving exact momentum routing, derivative kinematics, target provenance, and finite
equal-time regularisation. The current kinetic compiler supports only `gradient_order=Val(0)`.
""" wigner_transform

"""
    keldysh_component(x)

Return the Keldysh component of a propagator, self-energy, Fourier/Wigner result, or kinetic
self-energy returned by the public workflow.
"""
keldysh_component(x::DressedPropagator) = x.keldysh
keldysh_component(x::SelfEnergy) = x.keldysh
keldysh_component(x::FourierDressedPropagator) = x.keldysh
keldysh_component(x::FourierSelfEnergy) = x.keldysh
keldysh_component(x::WignerDressedPropagator) = x.keldysh
keldysh_component(x::WignerSelfEnergy) = x.keldysh
keldysh_component(x::KineticSelfEnergy) = x.keldysh

"""
    retarded_component(x)

Return the retarded component of a propagator, self-energy, Fourier/Wigner result, or kinetic
self-energy returned by the public workflow.
"""
retarded_component(x::DressedPropagator) = x.retarded
retarded_component(x::SelfEnergy) = x.retarded
retarded_component(x::FourierDressedPropagator) = x.retarded
retarded_component(x::FourierSelfEnergy) = x.retarded
retarded_component(x::WignerDressedPropagator) = x.retarded
retarded_component(x::WignerSelfEnergy) = x.retarded
retarded_component(x::KineticSelfEnergy) = x.retarded

"""
    advanced_component(x)

Return the advanced component of a propagator, self-energy, Fourier/Wigner result, or kinetic
self-energy returned by the public workflow.
"""
advanced_component(x::DressedPropagator) = x.advanced
advanced_component(x::SelfEnergy) = x.advanced
advanced_component(x::FourierDressedPropagator) = x.advanced
advanced_component(x::FourierSelfEnergy) = x.advanced
advanced_component(x::WignerDressedPropagator) = x.advanced
advanced_component(x::WignerSelfEnergy) = x.advanced
advanced_component(x::KineticSelfEnergy) = x.advanced

@doc """
    collision_offset(collision)

Return the distribution-independent part of a collision expression. For
`OffShellCollisionExpression` this is the exact full-frequency contribution
`I₀ = iΣᴷ`; no quasiparticle shell, occupation substitution, or finite-width prescription has
been applied.
""" collision_offset

@doc """
    collision_distribution_coefficient(collision)

Return the coefficient of the external statistical distribution in a collision expression. For
`OffShellCollisionExpression` this is the exact full-frequency contribution `I₁ = -A_Σ`, so
`I_coll = I₀ + F_target(k) I₁`. The external distribution is not absorbed into the internal-line
IR and the quasiparticle reduction remains a downstream operation.
""" collision_distribution_coefficient

@doc """
    kinematic_factor(sector)

Return the exact momentum polynomial multiplying a reduced collision sector.
""" kinematic_factor

@doc """
    frequency_support(sector)

Return the shell and principal-value frequency support carried by a reduced collision sector.
""" frequency_support

@doc """
    reduced_regular_terms(result)

Return the finite shell and principal-value terms after causal frequency reduction.
""" reduced_regular_terms

@doc """
    reduced_blocked_terms(result)

Return causal frequency structures that do not have a finite value in the strict reduction,
including genuine pinch singularities.
""" reduced_blocked_terms

@doc """
    reduced_trotter_terms(result)

Return finite-Trotter states that remain unresolved after structural equal-time reduction.
""" reduced_trotter_terms

@doc """
    collision_kernel_terms(kernel)

Return the canonical sector-to-occupation-polynomial mapping of a final `CollisionKernel`.
""" collision_kernel_terms

@doc """
    spectral_line_family(line)

Return the physical field family defining a finite-width spectral-line identity.
""" spectral_line_family

@doc """
    spectral_line_momentum(line)

Return the exact routed momentum defining a finite-width spectral-line identity.
""" spectral_line_momentum

@doc """
    spectral_data(model, line)

Return the explicit Lorentzian spectral data assigned to `line`. Missing line data is an error;
finite-width reduction never supplies an implicit cutoff or sharp-shell value.
""" spectral_data

@doc """
    spectral_energy(x[, line])

Inspect the dispersion center of Lorentzian spectral data, directly or through a line-resolved
spectral model.
""" spectral_energy

@doc """
    spectral_linewidth(x[, line])

Inspect the positive linewidth of Lorentzian spectral data, directly or through a line-resolved
spectral model.
""" spectral_linewidth

@doc """
    spectral_residue(x[, line])

Inspect the residue of Lorentzian spectral data, directly or through a line-resolved spectral
model.
""" spectral_residue

@doc """
    spectral_normalization(data)

Return the exact Lorentzian normalization `∫A dω/(2π) = Z`.
""" spectral_normalization

@doc """
    spectral_squared_weight(data)

Return the exact repeated-line integral `∫A² dω/(2π) = 2Z²/Γ`.
""" spectral_squared_weight

@doc """
    spectral_jacobian(weight_or_reduction)

Return the exact affine-frequency Jacobian of a resolved analytic Lorentzian spectral weight.
An unsupported reduction raises an `ArgumentError` rather than supplying a fictitious value.
""" spectral_jacobian

@doc """
    spectral_factors(weight_or_reduction)

Return the canonical repeated-line factors of a resolved analytic Lorentzian spectral weight.
An unsupported reduction raises an `ArgumentError`.
""" spectral_factors

@doc """
    spectral_reduction_kind(reduction)

Return `LorentzianResolved` or `LorentzianUnsupported` for a concrete finite-width reduction
result.
""" spectral_reduction_kind

@doc """
    spectral_reduction_resolved(reduction)

Return whether an analytic Lorentzian repeated-line reduction was resolved exactly.
""" spectral_reduction_resolved

@doc """
    convolution_reduction_kind(reduction)

Return `LorentzianConvolutionResolved` or `LorentzianConvolutionUnsupported` for a concrete
one-residual-shell reduction result.
""" convolution_reduction_kind

@doc """
    convolution_reduction_resolved(reduction)

Return whether the one-residual-shell Lorentzian convolution was resolved exactly.
""" convolution_reduction_resolved

@doc """
    convolution_jacobian(reduction)

Return the exact independent-frequency Jacobian of a resolved Lorentzian convolution.
""" convolution_jacobian

@doc """
    convolution_pivot_lines(reduction)

Return the spectral-line identities used as independent frequency coordinates.
""" convolution_pivot_lines

@doc """
    convolution_dependent_line(reduction)

Return the single residual spectral-line identity of a resolved Lorentzian convolution.
""" convolution_dependent_line

@doc """
    convolution_coefficients(reduction)

Return the exact rational coefficients expressing the dependent line frequency in the pivot
frequency coordinates.
""" convolution_coefficients

@doc """
    convolution_external_coefficient(reduction)

Return the exact coefficient multiplying the projected external frequency in the residual line.
""" convolution_external_coefficient

@doc """
    convolution_energy_mismatch(reduction)

Return the exact residual quasiparticle `EnergyForm`. Its zero-width limit is the strict
quasiparticle energy shell.
""" convolution_energy_mismatch

@doc """
    convolution_effective_linewidth(reduction, model)

Return `Γ_d + Σ |c_j|Γ_j` for a resolved one-residual-shell Lorentzian convolution.
""" convolution_effective_linewidth

@doc """
    spectral_line(x)

Return the physical spectral-line identity carried by a finite-width factor or power-counting
entry.
""" spectral_line

@doc """
    spectral_multiplicity(factor)

Return the repeated spectral multiplicity of one Lorentzian factor.
""" spectral_multiplicity

@doc """
    linewidth_exponent(x)

Return the linewidth exponent induced by one finite-width factor or power-counting entry.
""" linewidth_exponent

@doc """
    linewidth_powers(counting)

Return the explicit nonzero linewidth exponents retained separately from perturbative coupling
order in width-aware power counting.
""" linewidth_powers

@doc """
    width_aware_power_counting(parameter, weight_or_reduction)

Combine an exact perturbative `ParameterMonomial` with the linewidth powers implied by a
resolved Lorentzian spectral weight. The linewidth scaling is not collapsed into coupling order.
""" width_aware_power_counting

@doc """
    lorentzian_spectral_data_from_retarded_self_energy(ε, Σᴿ; residue=1)

Construct Lorentzian line data from an on-shell retarded self-energy with
`E = ε + Re Σᴿ` and `Γ = -2 Im Σᴿ`. The quasiparticle residue is supplied explicitly and defaults
to one.
""" lorentzian_spectral_data_from_retarded_self_energy

@doc """
    lorentzian_spectral_value(data, ω)

Evaluate the explicit Lorentzian spectral density.
""" lorentzian_spectral_value

@doc """
    lorentzian_integrated_power(data, multiplicity)

Evaluate the exact integrated positive integer power of a Lorentzian spectral line.
""" lorentzian_integrated_power

@doc """
    lorentzian_spectral_reduction(term)

Attempt the exact analytic Lorentzian frequency reduction for a factorized full-rank
all-spectral `SpectralDispersiveTerm`. The function always returns a concrete
`LorentzianSpectralReduction`; unsupported dispersive, shifted, rank-deficient, or nonfactorized
structures are tagged `LorentzianUnsupported` and are not regularised implicitly.
""" lorentzian_spectral_reduction

@doc """
    lorentzian_convolution_reduction(term, target)

Attempt the exact analytic `L + 1` line Lorentzian convolution over `L` loop frequencies. The
result is a concrete `LorentzianConvolutionReduction` and never introduces numerical quadrature.
""" lorentzian_convolution_reduction

@doc """
    evaluate_spectral_weight(weight_or_reduction, model)

Evaluate a resolved factorized Lorentzian spectral weight using explicit line-resolved data.
Every required line must be present in the model, and unsupported reductions are rejected.
""" evaluate_spectral_weight

@doc """
    evaluate_spectral_convolution(reduction, model, ω_external)

Evaluate a resolved one-residual-shell Lorentzian convolution at the supplied external energy.
Every internal line must have explicit spectral data.
""" evaluate_spectral_convolution
