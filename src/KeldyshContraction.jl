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
    include("propagator_algebra/canonicalize.jl")
    include("propagator_algebra/diagram.jl")
    include("dressed_propagator.jl")

    include("wick_contractions.jl")
    include("filters.jl")
    include("self_energy.jl")

    # Statistics extensions
    include("fermionic_keldysh.jl")
    include("fermionic_propagator_dispatch.jl")

    include("propagator_algebra/wigner.jl")
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
        wigner_transform,
        parameters
end

end
