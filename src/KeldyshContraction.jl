"""
$(DocStringExtensions.README)
"""
module KeldyshContraction

using DocStringExtensions: DocStringExtensions
using LinearAlgebra: LinearAlgebra
using EnumX: @enumx

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

    # Legacy private routing used only by the pre-Fourier collision reducer.
    include("propagator_algebra/legacy_momentum_routing.jl")
    include("collision_integral.jl")

    # show methods
    include("show_methods/latexify_recipes.jl")
    include("show_methods/printing.jl")

    export @qfields,
        FieldFamily,
        Field,
        field_family,
        field_families,
        target_family,
        ParameterMonomial,
        parameter_monomial,
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
        FourierDiagram,
        FourierDiagrams,
        FourierDressedPropagator,
        FourierSelfEnergy,
        wigner_transform,
        WignerDiagram,
        WignerDiagrams,
        WignerDressedPropagator,
        WignerSelfEnergy,
        HomogeneousWignerContext,
        gradient_order,
        external_wigner_momentum,
        wigner_context,
        KineticLineKind,
        KineticRetarded,
        KineticAdvanced,
        KineticSpectral,
        StatisticalWeight,
        NoStatisticalWeight,
        DistributionWeight,
        KineticLine,
        KineticMonomial,
        KineticTerm,
        KineticExpression,
        KineticSelfEnergy,
        kinetic_expression,
        kinetic_lines,
        kinetic_line_kind,
        statistical_weight,
        regularisation_shift,
        retarded_minus_advanced,
        spectral_self_energy,
        statistical_occupation_coefficients,
        statistical_from_occupation,
        OffShellCollisionExpression,
        off_shell_collision_expression,
        collision_offset,
        collision_distribution_coefficient,
        parameters
end

end
