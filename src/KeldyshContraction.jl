"""
$(DocStringExtensions.README)
"""
module KeldyshContraction

using DocStringExtensions

using TermInterface: TermInterface
using SymbolicUtils: SymbolicUtils, @syms, arguments
using Combinatorics: Combinatorics
using SmallCollections: SmallCollections
using StaticArrays: StaticArrays, SVector

using Latexify
using MacroTools: MacroTools
using LaTeXStrings

# utils
include("utils.jl")

# Fields
include("keldysh_algebra/interface.jl")
include("keldysh_algebra/keldysh_algebra.jl")
include("keldysh_algebra/QTerm.jl")
include("keldysh_algebra/field_math.jl")
include("keldysh_algebra/hashing.jl")
include("InteractionLagrangian.jl")

# Propagators
include("propagator.jl")
include("diagram.jl")
include("dressed_propagator.jl")

include("wick_contractions.jl")
include("filters.jl")
include("canonicalize.jl")
include("self_energy.jl")

# show methods
include("show_methods/latexify_recipes.jl")
include("show_methods/printing.jl")

export @qfields,
    Destroy,
    Create,
    wick_contraction,
    Quantum,
    Classical,
    DressedPropagator,
    SelfEnergy,
    InteractionLagrangian,
    @syms,
    arguments,
    topologies

end
