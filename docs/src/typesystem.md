# Making symbolic quantum field algebra

The implemented Computer Algebra System (CAS) for the quantum fields is build using [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl) and [TermInterface.jl](https://github.com/JuliaSymbolics/TermInterface.jl/):

- [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl) is a utility library for symbolic computation in Julia. It provides a set of tools for manipulating and transforming symbolic expressions, including support for rule-base algebraic operations, simplification, and pattern matching.
- [TermInterface.jl](https://github.com/JuliaSymbolics/TermInterface.jl/) is the package that defines the common interface for defining symbolic expressions. It allows users to define their own types and "define" them as symbols of an algebra.

Using these packages, we can define the algebra of quantum fields in a symbolic way. For this we proposed the type hierarchy:

```@example interface
using GraphRecipes, Plots, KeldyshContraction, Random
Random.seed!(1) # hide
theme(:dracula) # hide
plot(
  KeldyshContraction.QField;
  method=:tree, fontsize=10, markersize = 0.12, nodeshape=:ellipse
)
```

[`KeldyshContraction.QSym`](@ref) will be abstract type representing the individual field of type [`Destroy`](@ref) and [`Create`](@ref). [`KeldyshContraction.QTerm`](@ref) will represent the terms of the algebra, which are the *products* and *sum* of the fields. The type naming and hierarchy is heavily inspired by the implementation in `QuantumCumulants.jl`.

QSym can then have additional properties to make it Keldsysh-specific type fields. In the package this is done by adding fields to `Create` and `Destroys` using *Enum* objects:

- [`KeldyshContraction.KeldyshContour`](@ref) - the Keldysh contour of the field, which can be either `KeldyshContour.Quantum` or `KeldyshContour.Classical`.
- [`KeldyshContraction.Regularisation`](@ref) - the tadpole regularisation of the field, which can be either `Regularisation.Zero`, `Regularisation.Plus` or `Regularisation.Minus`.
- [`KeldyshContraction.Position`](@ref) - the position of the field, which can be either `Position.In`, `Position.Out` or `Position.Bulk`.

To make our quantum field types work with the symbolic algebra system, we need to implement several interface functions from SymbolicUtils.jl and TermInterface.jl:

```julia
using KeldyshContraction: SymbolicUtils, TermInterface

@qfields ϕ::Destroy(Classical)

TermInterface.head(ϕ) = :call
SymbolicUtils.iscall(ϕ) = false
SymbolicUtils.iscall(ϕ*ϕ) = true
TermInterface.metadata(ϕ) = nothing
SymbolicUtils.symtype(ϕ) = Destroy{Classical, KeldyshContraction.Zero, Nothing} 
(one(ϕ), zero(ϕ)) = (1, 0)
```

The key interfaces are:

- `head`: Defines how the expression should be interpreted (`:call` indicates function application)
- `iscall`: Specifies which types represent function calls (`QTerm` types) vs atomic symbols (`QSym` types)
- `metadata`: Allows attaching additional information to terms in this package is not used and set to nothing

For type promotion during operations, we implement:

```julia
# Type promotion rules for operations like +, -, *, /, //, \, ^
SymbolicUtils.promote_symtype(::typeof(+), T::Type{<:QField}, S::Type{<:QField}) = promote_type(T, S)
SymbolicUtils.promote_symtype(::typeof(*), T::Type{<:QField}, S::Type{<:Number}) = T
```

These implementations allow our quantum field types to:

1. Be recognized as symbolic terms
2. Participate in algebraic operations
3. Follow proper type promotion rules
4. Handle basic mathematical concepts like identities and zeros
