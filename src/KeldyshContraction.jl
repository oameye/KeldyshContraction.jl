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

using Latexify: Latexify, latexify, @latexrecipe
using MacroTools: MacroTools
using LaTeXStrings: LaTeXStrings, latexstring

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
include("propagator_algebra/propagator.jl")
include("propagator_algebra/canonicalize.jl")
include("propagator_algebra/Momentum.jl")
include("propagator_algebra/diagram.jl")
include("propagator_algebra/wigner.jl")

include("dressed_propagator.jl")

include("wick_contractions.jl")
include("filters.jl")
include("self_energy.jl")

# show methods
include("show_methods/latexify_recipes.jl")
include("show_methods/printing.jl")

const Quantum = KeldyshContour.Quantum
const Classical = KeldyshContour.Classical

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
