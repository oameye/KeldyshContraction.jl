"""
    KeldyshContraction

Symbolic Schwinger--Keldysh perturbation theory for bosonic and fermionic fields, from
interaction actions and 1PI self-energies to quantum-kinetic collision kernels.
"""
module KeldyshContraction

using DocStringExtensions: DocStringExtensions
using LinearAlgebra: LinearAlgebra
using EnumX: @enumx
using SciMLPublic: @public

using TermInterface: TermInterface
using SymbolicUtils: SymbolicUtils, @syms, arguments
using Combinatorics: Combinatorics
using SmallCollections: SmallCollections, FixedVector
using NautyGraphs: NautyGraphs
using Graphs: Graphs
using DispatchDoctor: @stable

const ComplexRationals = Complex{Rational{Int64}}

@stable default_mode = "disable" default_codegen_level = "min" begin
    # utils
    include("utils.jl")

    # Fields
    include("keldysh_algebra/interface.jl")
    include("keldysh_algebra/derivative_metadata.jl")
    include("keldysh_algebra/keldysh_algebra.jl")
    include("keldysh_algebra/QTerm.jl")
    include("keldysh_algebra/field_math.jl")
    include("keldysh_algebra/hashing.jl")
    include("parameters.jl")
    include("InteractionLagrangian.jl")

    # Propagators
    include("propagator_algebra/propagator.jl")
    include("propagator_algebra/momentum_algebra.jl")
    include("propagator_algebra/momentum_routing.jl")
    include("propagator_algebra/momentum_affine_routing.jl")
    include("propagator_algebra/canonicalize.jl")
    include("propagator_algebra/diagram.jl")
    include("equal_time_regularisation.jl")
    include("propagator_algebra/fourier_diagram.jl")
    include("propagator_algebra/momentum_polynomial.jl")
    include("dressed_propagator.jl")

    include("wick_contractions.jl")
    include("filters.jl")
    include("self_energy.jl")

    # Statistics extensions
    include("fermionic_keldysh.jl")
    include("fermionic_propagator_dispatch.jl")

    include("propagator_algebra/fourier_transform.jl")
    include("propagator_algebra/wigner.jl")
    include("propagator_algebra/spectral_statistical.jl")
    include("off_shell_collision.jl")
    include("spectral_dispersive_collision.jl")
    include("frequency_support.jl")
    include("full_rank_frequency_reduction.jl")
    include("finite_width_spectral.jl")
    include("finite_width_convolution.jl")
    include("finite_width_spectral_variation.jl")
    include("finite_width_frequency_collision.jl")
    include("finite_width_collision_evaluation.jl")
    include("canonical_frequency_expression.jl")
    include("causal_frequency_integration.jl")
    include("affine_singular_support.jl")
    include("causal_frequency_support_reduction.jl")
    include("exceptional_frequency_classification.jl")
    include("statistical_frequency_algebra.jl")
    include("collision_statistical_algebra.jl")
    include("occupation_background.jl")
    include("occupation_linearization.jl")
    include("canonical_frequency_collision.jl")
    include("trotter_frequency_reduction.jl")
    include("trotter_frequency_collision.jl")
    include("collision_frequency_assembly.jl")
    include("spectral_self_energy_kernel.jl")
    include("spectral_self_energy_background.jl")
    include("finite_width_occupation_reduction.jl")
    include("finite_width_collision_response.jl")
    include("collision_projection_core.jl")
    include("finite_width_projected_matrix.jl")
    include("collision_momentum_quotient.jl")
    include("collision_momentum_projective.jl")
    include("collision_kernel_api.jl")

    # Legacy private routing used only by the pre-Fourier collision reducer.
    include("propagator_algebra/legacy_momentum_routing.jl")
    include("collision_integral.jl")

    # show methods
    include("show_methods/latexify_recipes.jl")
    include("show_methods/printing.jl")
    include("show_methods/kinetic_pipeline.jl")
    include("show_methods/physical_pipeline.jl")
    include("show_methods/physical_compaction.jl")
    include("show_methods/kernel_display_gauge.jl")

    # Stable public inspection helpers defined after the full representation stack exists.
    include("public_api_accessors.jl")

    # Qualified public API is deliberately small: semantic inspection and the stable physical
    # boundaries of the kinetic compiler. Lower-level representation machinery remains private.
    @public Regularisation,
    OffShellCollisionExpression,
    SpectralLineIdentity,
    LorentzianSpectralData,
    LorentzianSpectralModel,
    LorentzianSpectralDataVariation,
    LorentzianSpectralModelVariation,
    LorentzianSpectralFactor,
    LorentzianSpectralWeight,
    LorentzianSpectralReductionKind,
    LorentzianSpectralReduction,
    LorentzianResolved,
    LorentzianUnsupported,
    LorentzianConvolutionReductionKind,
    LorentzianConvolutionReduction,
    LorentzianConvolutionResolved,
    LorentzianConvolutionUnsupported,
    LinewidthPower,
    WidthAwarePowerCounting,
    FiniteWidthFrequencyReductionKind,
    FiniteWidthFactorized,
    FiniteWidthConvolution,
    FiniteWidthUnsupported,
    FiniteWidthFrequencyTerm,
    FiniteWidthFrequencyExpression,
    FiniteWidthFrequencyCollision,
    FiniteWidthEvaluatedExpression,
    FiniteWidthEvaluatedCollision,
    FiniteWidthOccupationCollision,
    OccupationBackground,
    OccupationLinearization,
    BackgroundOccupationLinearization,
    OccupationSpectralResponse,
    FiniteWidthCollisionLinearization,
    SpectralSelfEnergySector,
    SpectralSelfEnergyKernel,
    BackgroundSpectralSelfEnergyKernel,
    CollisionProjectionBasis,
    CollisionPerturbationClosure,
    CollisionProjectionFunctional,
    ProjectedCollisionMatrix,
    field_families,
    target_family,
    parameters,
    matrix,
    order,
    statistics,
    gradient_order,
    wigner_context,
    keldysh_component,
    retarded_component,
    advanced_component,
    collision_offset,
    collision_distribution_coefficient,
    momentum_basis,
    external_wigner_momentum,
    kinematic_factor,
    frequency_support,
    reduced_regular_terms,
    reduced_blocked_terms,
    reduced_trotter_terms,
    collision_kernel_terms,
    spectral_self_energy_terms,
    background_spectral_self_energy_terms,
    spectral_self_energy_blocked_terms,
    spectral_self_energy_trotter_terms,
    evaluate_spectral_self_energy_background,
    spectral_line_family,
    spectral_line_momentum,
    spectral_data,
    spectral_energy,
    spectral_linewidth,
    spectral_residue,
    spectral_data_variation,
    spectral_energy_variation,
    spectral_linewidth_variation,
    spectral_residue_variation,
    spectral_normalization,
    spectral_squared_weight,
    spectral_jacobian,
    spectral_factors,
    spectral_reduction_kind,
    spectral_reduction_resolved,
    convolution_reduction_kind,
    convolution_reduction_resolved,
    convolution_jacobian,
    convolution_pivot_lines,
    convolution_dependent_line,
    convolution_coefficients,
    convolution_external_coefficient,
    convolution_energy_mismatch,
    convolution_effective_linewidth,
    spectral_line,
    spectral_multiplicity,
    linewidth_exponent,
    linewidth_powers,
    width_aware_power_counting,
    finite_width_source_term,
    finite_width_reduction_kind,
    finite_width_spectral_reduction,
    finite_width_convolution_reduction,
    finite_width_frequency_collision,
    finite_width_resolved_terms,
    finite_width_unsupported_terms,
    evaluate_finite_width_collision,
    finite_width_occupation_terms,
    finite_width_unsupported_offset_terms,
    finite_width_unsupported_distribution_terms,
    finite_width_occupation_collision,
    occupation_background_value,
    evaluate_occupation_polynomial,
    occupation_linearization_terms,
    occupation_linearization,
    background_occupation_linearization_terms,
    evaluate_occupation_linearization,
    spectral_occupation_response_terms,
    finite_width_linearized_terms,
    finite_width_occupation_response_terms,
    finite_width_spectral_response_terms,
    linearize_finite_width_collision,
    left_projection_basis,
    right_perturbation_basis,
    perturbation_amplitude,
    project_collision_channel,
    projected_collision_matrix,
    lorentzian_spectral_data_from_retarded_self_energy,
    lorentzian_spectral_value,
    lorentzian_integrated_power,
    lorentzian_integrated_power_variation,
    lorentzian_spectral_reduction,
    lorentzian_convolution_reduction,
    evaluate_spectral_weight,
    evaluate_spectral_convolution,
    evaluate_spectral_weight_variation,
    evaluate_spectral_convolution_variation,
    spectral_self_energy_kernel

    # Small workflow surface intended for ordinary unqualified use.
    export @qfields,
        FieldFamily,
        Field,
        field_family,
        Boson,
        Fermion,
        bar,
        partial,
        derivatives,
        wick_contraction,
        Quantum,
        Classical,
        One,
        Two,
        DressedPropagator,
        SelfEnergy,
        InteractionLagrangian,
        convert_coefficients,
        rationalize_coefficients,
        @syms,
        arguments,
        topologies,
        fourier_transform,
        wigner_transform,
        kinetic_expression,
        off_shell_collision_expression,
        spectral_dispersive_collision,
        reduce_frequency_collision,
        occupation_reduced_expression,
        quotient_loop_momenta,
        CollisionKernel,
        collision_kernel
end

end
