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

The package uses one concrete field family, parameterized only by statistics. Bosonic
fields therefore have type `Field{Boson}`; barred/unbarred orientation, Keldysh component,
position, regularisation, and internal indices are concrete value data.

```@docs
KeldyshContraction.QField
KeldyshContraction.QSym
KeldyshContraction.Field
KeldyshContraction.Boson
KeldyshContraction.bar
```

#### Field properties

Bosonic `Classical` and `Quantum` are semantic labels over the package's neutral two-valued
Keldysh index. They do not create different Julia field types.

```@docs
KeldyshContraction.Regularisation
KeldyshContraction.Position
KeldyshContraction.In
KeldyshContraction.Out
KeldyshContraction.Bulk
```

#### Field constructors

Fields are normally created with `@qfields`:

```@example API
using KeldyshContraction
using KeldyshContraction: position

@qfields ϕ::Boson(Classical)
barϕ = bar(ϕ)

(position(ϕ), typeof(ϕ), typeof(barϕ))
```

```@docs
KeldyshContraction.@qfields
```

Calling a field with a `Position` or `Regularisation` value returns the same concrete field
type with that value changed.

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
arguments(::KeldyshContraction.QMul)
arguments(::KeldyshContraction.QAdd)
KeldyshContraction.convert_coefficients
KeldyshContraction.rationalize_coefficients
```

Coefficient conversion is explicit. In particular, constructing an
`InteractionLagrangian` does not rationalize floating-point coefficients according to their
runtime values.

The properties of an expression can be checked using:

```@docs
KeldyshContraction.is_bulk
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

```@docs
KeldyshContraction.PropagatorType
DressedPropagator
KeldyshContraction.matrix(::DressedPropagator)
```

### Self-energy

```@docs
KeldyshContraction.SelfEnergy
KeldyshContraction.matrix(::KeldyshContraction.SelfEnergy)
```
