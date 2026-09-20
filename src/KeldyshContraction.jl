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
import GraphCombinations as GC
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
    include("gc_wick_backend.jl")
    include("self_energy.jl")

    # Statistics extensions
    include("fermionic_keldysh.jl")
    include("gc_wick_causal_pruning.jl")
    include("gc_self_energy_pruning.jl")
    include("gc_wick_routing.jl")
    include("fermionic_propagator_dispatch.jl")

    include("propagator_algebra/fourier_transform.jl")
    include("propagator_algebra/wigner.jl")
    include("propagator_algebra/spectral_statistical.jl")
    include("off_shell_collision.jl")
    include("spectral_dispersive_collision.jl")
    include("frequency_support.jl")
    include("full_rank_frequency_reduction.jl")
    include("canonical_frequency_expression.jl")
    include("causal_frequency_integration.jl")
    include("affine_singular_support.jl")
    include("causal_frequency_support_reduction.jl")
    include("exceptional_frequency_classification.jl")
    include("statistical_frequency_algebra.jl")
    include("collision_statistical_algebra.jl")
    include("canonical_frequency_collision.jl")
    include("trotter_frequency_reduction.jl")
    include("trotter_frequency_collision.jl")
    include("collision_frequency_assembly.jl")
    include("collision_momentum_quotient.jl")
    include("collision_momentum_projective.jl")
    include("collision_momentum_graphcombinations.jl")
    include("collision_momentum_workspace.jl")
    include("collision_kernel_api.jl")

    # Legacy private routing used only by the pre-Fourier collision reducer.
    include("propagator_algebra/legacy_momentum_routing.jl")
    include("collision_integral.jl")

    # show methods
    include("show_methods/latexify_recipes.jl")
    include("show_methods/printing.jl")
    include("show_methods/kinetic_pipeline.jl")

    # Stable public inspection helpers defined after the full representation stack exists.
    include("public_api_accessors.jl")

    # Qualified public API is deliberately small: only semantic inspection that users need
    # in the manual or canonical examples. Representation and compiler IR remain private.
    @public Regularisation,
    field_families,
    target_family,
    parameters,
    matrix,
    order,
    statistics,
    keldysh_component,
    retarded_component,
    advanced_component,
    kinematic_factor,
    frequency_support,
    reduced_regular_terms,
    reduced_blocked_terms,
    reduced_trotter_terms,
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
        LoopMomentumQuotientWorkspace,
        quotient_loop_momenta,
        CollisionKernel,
        collision_kernel
end

end