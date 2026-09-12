"""
$(DocStringExtensions.README)
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
    include("causal_frequency_reduction.jl")
    include("higher_order_causal_frequency_reduction.jl")
    include("full_rank_frequency_reduction.jl")
    include("dependent_shell_support.jl")
    include("partial_spectral_frequency_reduction.jl")
    include("exceptional_frequency_classification.jl")
    include("trotter_frequency_reduction.jl")
    include("state_frequency_reduction.jl")
    include("state_exceptional_frequency.jl")
    include("structured_causal_frequency_reduction.jl")
    include("reduced_frequency_result.jl")
    include("collision_statistical_algebra.jl")
    include("collision_frequency_assembly.jl")
    include("collision_momentum_quotient.jl")
    include("collision_momentum_projective.jl")

    # Legacy private routing used only by the pre-Fourier collision reducer.
    include("propagator_algebra/legacy_momentum_routing.jl")
    include("collision_integral.jl")

    # show methods
    include("show_methods/latexify_recipes.jl")
    include("show_methods/printing.jl")

    # Stable qualified API. These names are intentionally public without being imported by
    # `using KeldyshContraction`; SciMLPublic backports Julia's `public` declaration to 1.10.
    @public Statistics,
        QField,
        QSym,
        KeldyshIndex,
        Orientation,
        Regularisation,
        Position,
        IndexKind,
        Bulk,
        In,
        Out,
        reconstruct,
        position,
        field_families,
        target_family,
        ParameterMonomial,
        parameter_monomial,
        parameters,
        QTerm,
        QMul,
        QAdd,
        coefficient,
        fields,
        terms,
        exchange_sign,
        is_conserved,
        is_physical,
        LagrangianSum,
        DressedPropagatorSum,
        PropagatorType,
        matrix,
        SelfEnergySum,
        MomentumVariable,
        LinearMomentum,
        MomentumBasis,
        MomentumComponent,
        MomentumMonomial,
        MomentumPolynomial,
        momentum_basis,
        kinematic_factor,
        FourierDiagram,
        FourierDiagrams,
        FourierDressedPropagator,
        FourierSelfEnergy,
        WignerDiagram,
        WignerDiagrams,
        WignerDressedPropagator,
        WignerSelfEnergy,
        HomogeneousWignerContext,
        gradient_order,
        external_wigner_momentum,
        wigner_context,
        KineticSelfEnergy,
        retarded_minus_advanced,
        spectral_self_energy,
        statistical_occupation_coefficients,
        statistical_from_occupation,
        OffShellCollisionExpression,
        collision_offset,
        collision_distribution_coefficient,
        SpectralDispersiveCollision,
        EnergyShell,
        PrincipalValueSupport,
        FrequencySupport,
        frequency_support,
        frequency_factor,
        ReducedFrequencyCollision,
        reduced_regular_terms,
        reduced_dependent_terms,
        reduced_causal_terms,
        reduced_trotter_terms,
        OccupationReducedExpression,
        occupation_reduced_terms,
        LoopQuotientedExpression,
        loop_quotient_terms,
        CollisionKernelSector,
        collision_kernel_terms

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
