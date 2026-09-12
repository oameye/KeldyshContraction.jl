```@meta
CollapsedDocStrings = true
```

# API

```@contents
Pages = ["API.md"]
Depth = 2:3
```

```@docs
KeldyshContraction
```

## Field types

### Individual fields

```@example API
using Term, KeldyshContraction # hide
Term.typestree(KeldyshContraction.QSym) # hide
```

The package uses one concrete field representation, parameterized only by statistics.
Bosonic and fermionic fields therefore have types `Field{Boson}` and `Field{Fermion}` while
sharing the same storage layout. Physical field identity is stored in `FieldFamily{S}`;
barred/unbarred orientation, Keldysh component, position, regularisation, and coordinate-
derivative decoration are concrete value data on each field component.

```@docs
KeldyshContraction.QField
KeldyshContraction.QSym
KeldyshContraction.Statistics
KeldyshContraction.Boson
KeldyshContraction.Fermion
KeldyshContraction.FieldFamily
KeldyshContraction.Field
KeldyshContraction.field_family
KeldyshContraction.bar
KeldyshContraction.partial
KeldyshContraction.derivatives
```

#### Field properties

The stored Keldysh index is neutral and two-valued. Bosonic `Classical` / `Quantum` and
fermionic `One` / `Two` are statistics-specific semantic labels over that same representation;
none of these component labels creates a different Julia field type.

For fermions, `bar(psi)` constructs the independent barred Grassmann path-integral variable
and preserves the `One` / `Two` component. It is not an operator-adjoint operation. See
[the convention page](conventions.md) for the asymmetric fermionic LO rotation and matrix
placement.

Coordinate derivatives likewise preserve the outer field type:

```@example API
using KeldyshContraction

@qfields ϕd::Boson ψd::Fermion
c = ϕd[Classical]
ψ1 = ψd[One]

∂xc = partial(c, :x)
∂xyψ = partial(partial(ψ1, :y), :x)

(typeof(∂xc), typeof(∂xyψ), derivatives(∂xyψ))
```

The supported coordinate axes are `:t`, `:x`, `:y`, and `:z`; unsupported axis labels raise
`ArgumentError`. Derivative axes are stored canonically as concrete value-level metadata.
Coordinate derivatives commute as operations on the field label, so
`partial(partial(ψ1,:x),:y)` and `partial(partial(ψ1,:y),:x)` denote the same generator.
`derivatives(field)` returns an owned copy of the canonical axis sequence.

Derivative decoration is part of generator and diagram identity but not of the underlying
`FieldFamily`. Consequently `ψ` and `∂xψ` are distinct Grassmann generators of the same
physical field family. This coordinate-space layer deliberately does not convert derivatives
to momentum factors; that belongs to the Fourier/momentum transformation.

```@docs
KeldyshContraction.KeldyshIndex
KeldyshContraction.Orientation
KeldyshContraction.Regularisation
KeldyshContraction.Position
KeldyshContraction.IndexKind
KeldyshContraction.Bulk
KeldyshContraction.In
KeldyshContraction.Out
KeldyshContraction.reconstruct
```

The constructors `Bulk(i)`, `In()`, and `Out()` create `Position` values. Calling a field
with a `Position` or `Regularisation` value returns the same concrete field type with that
value changed and preserves its field family and derivative decoration.

#### Field constructors

Declare physical field families with `@qfields`, then obtain their Keldysh components by
indexing the family:

```@example API
using KeldyshContraction
using KeldyshContraction: position

@qfields ϕ::Boson ψ::Fermion
c, q = ϕ[Classical], ϕ[Quantum]
ψ1, ψ2 = ψ[One], ψ[Two]
barψ1 = bar(ψ1)

(
    field_family(c) == field_family(q),
    field_family(ψ1) == field_family(ψ2),
    position(c),
    typeof(c),
    typeof(ψ1),
    typeof(barψ1),
)
```

```@docs
KeldyshContraction.@qfields
```

### Field algebra

```@example API
using Term, KeldyshContraction # hide
Term.typestree(KeldyshContraction.QTerm) # hide
```

Products and sums are homogeneous concrete containers:

```text
QMul{C,S} -> Vector{Field{S}}
QAdd{C,S} -> Vector{QMul{C,S}}
```

where `C` is the coefficient representation and `S` the field statistics. Algebraic zero
and one remain inside this symbolic representation instead of returning value-dependent raw
scalars. For `Fermion`, canonical products additionally implement Grassmann exchange signs and
nilpotency while retaining the same `QMul{C,Fermion}` result representation.

Derivative decoration participates in that same algebra. In particular, for a fermion field
`ψ`, `ψ*ψ` and `partial(ψ,:x)*partial(ψ,:x)` vanish, while `ψ*partial(ψ,:x)` is generally
nonzero and changes sign under exchange. Bosonic differentiated generators retain bosonic
commutation.

```@docs
KeldyshContraction.QTerm
KeldyshContraction.QMul
KeldyshContraction.QAdd
KeldyshContraction.coefficient
KeldyshContraction.fields
KeldyshContraction.terms
KeldyshContraction.convert_coefficients
KeldyshContraction.rationalize_coefficients
KeldyshContraction.exchange_sign
```

`SymbolicUtils.arguments` remains available for symbolic-tree interoperability. For package
code, use the semantic accessors `coefficient`, `fields`, and `terms` instead of depending
on the mixed SymbolicUtils argument vector.

Coefficient conversion is explicit. In particular, constructing an
`InteractionLagrangian` does not rationalize floating-point coefficients according to their
runtime values.

```@docs
KeldyshContraction.is_conserved
KeldyshContraction.is_physical
```

### Perturbation parameters

Perturbative result dictionaries use a package-native canonical commutative monomial as
their key. SymbolicUtils symbols and products are accepted at API boundaries and normalized
immediately, so symbolic expression trees do not become part of the computational result
types.

```@docs
KeldyshContraction.ParameterMonomial
KeldyshContraction.parameter_monomial
```

## Systems

```@docs
InteractionLagrangian
KeldyshContraction.LagrangianSum
```

Single-family interactions infer their target propagator. Multi-family interactions use the
same explicit `target` field-family selection for bosonic and fermionic statistics.

Adding `InteractionLagrangian`s with common physical field families constructs a
`LagrangianSum{C,S}`. `DressedPropagator` supports these sums for both `Boson` and `Fermion`;
distinct perturbation processes remain separated by canonical `ParameterMonomial` keys, and
higher orders include the corresponding mixed monomials.

## Wick contraction

The perturbation order and propagator edge count are supplied as `Val` arguments because
they determine the static diagram representation. Fermionic contractions reuse the same
`WickPairing{S,E}` representation; statistics dispatch supplies permutation parity to the
pairing weight. Coordinate derivatives are passive for pairing parity and R/A/K
classification but remain attached to their exact contraction endpoints.

```@docs
wick_contraction
```

The coordinate-space result representation is stable but intentionally qualified:
`KeldyshContraction.Diagrams` is an iterable collection of `KeldyshContraction.Diagram`
objects, whose propagator edges are `KeldyshContraction.Edge` values. Use the public
`contractions`, `topology`, `fields`, and `propagator_type` accessors for inspection rather
than depending on storage fields. `topologies` remains the exported grouping helper.

### Propagator

Propagator edges carry their retarded, advanced, Keldysh, or spectral component as concrete
value data. `DressedPropagator` stores semantic R/A/K components independently of statistics;
`matrix` supplies the statistics-specific layout. In particular, the fermionic LO matrix is
`[[R,K],[0,A]]`, while the bosonic RAK matrix is `[[K,R],[A,0]]`.

For a `LagrangianSum`, `DressedPropagator` returns a parameter-keyed
`DressedPropagatorSum`. Each stored value is the same concrete `DressedPropagator{C,S,...}`
used by the single-interaction path. Use `keldysh_component`, `retarded_component`, and
`advanced_component` rather than depending on the storage fields. The qualified `order`,
`statistics`, `parameters`, and `target_family` accessors expose provenance without widening
the default namespace.

```@docs
DressedPropagator
KeldyshContraction.DressedPropagatorSum
KeldyshContraction.PropagatorType
KeldyshContraction.matrix(::DressedPropagator)
KeldyshContraction.keldysh_component
KeldyshContraction.retarded_component
KeldyshContraction.advanced_component
```

### Self-energy

`SelfEnergy` likewise stores semantic R/A/K components. The fermionic LO self-energy matrix
is `[[R,K],[0,A]]`; the bosonic self-energy uses `[[0,A],[R,K]]`.

Applying `SelfEnergy` to a `DressedPropagatorSum` preserves the parameter-monomial keys and
returns a concrete `SelfEnergySum`. Derivative endpoint decoration is preserved through
self-energy extraction and remains coordinate-space metadata until the later Fourier layer
converts it to a momentum polynomial. The same component and provenance accessors used by
`DressedPropagator` apply to `SelfEnergy`.

```@docs
KeldyshContraction.SelfEnergy
KeldyshContraction.SelfEnergySum
KeldyshContraction.matrix(::KeldyshContraction.SelfEnergy)
```
