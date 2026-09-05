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
Bosonic fields therefore have type `Field{Boson}`. Physical field identity is stored in a
`FieldFamily{Boson}`, while barred/unbarred orientation, Keldysh component, position, and
regularisation are concrete value data on each field component.

```@docs
KeldyshContraction.QField
KeldyshContraction.QSym
KeldyshContraction.Statistics
KeldyshContraction.Boson
KeldyshContraction.FieldFamily
KeldyshContraction.Field
KeldyshContraction.field_family
KeldyshContraction.bar
```

#### Field properties

Bosonic `Classical` and `Quantum` are semantic aliases over the package's neutral two-valued
Keldysh index. They are components of one physical field family and do not create different
Julia field types.

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
value changed and preserves its field family.

#### Field constructors

Declare physical field families with `@qfields`, then obtain their Keldysh components by
indexing the family:

```@example API
using KeldyshContraction
using KeldyshContraction: position

@qfields ϕ::Boson
c, q = ϕ[Classical], ϕ[Quantum]
barc = bar(c)

(field_family(c) == field_family(q), position(c), typeof(c), typeof(barc))
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
and one remain inside this symbolic representation instead of returning value-dependent
raw scalars.

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

## Systems

```@docs
InteractionLagrangian
```

## Wick contraction

The perturbation order and propagator edge count are supplied as `Val` arguments because
they determine the static diagram representation.

```@docs
wick_contraction
```

### Propagator

Propagator edges carry their retarded, advanced, Keldysh, or spectral component as concrete
value data.

```@docs
DressedPropagator
KeldyshContraction.PropagatorType
KeldyshContraction.matrix(::DressedPropagator)
```

### Self-energy

```@docs
KeldyshContraction.SelfEnergy
KeldyshContraction.matrix(::KeldyshContraction.SelfEnergy)
```
