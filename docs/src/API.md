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

## Field Types

### Individual Fields

```@example API
using Term, KeldyshContraction # hide
Term.typestree(KeldyshContraction.QSym) # hide
```

```@docs
KeldyshContraction.QField
KeldyshContraction.QSym
KeldyshContraction.Create
KeldyshContraction.Destroy
```

#### Field Properties

The field properties are determined by the Enum objects:

```@docs
KeldyshContraction.KeldyshContour
KeldyshContraction.Regularisation
```

And the position of the field is determined by the `AbstractPosition` object:

```@docs
KeldyshContraction.Position
KeldyshContraction.In
KeldyshContraction.Out
KeldyshContraction.Bulk
```

#### Field Constructors

It is expected to create the fields using the `Create` and `Destroy` constructors, together with the `@qfields` macro.

```@docs
KeldyshContraction.@qfields
```

The created fields are callable to change a property of the individual fields:

```@example API
using KeldyshContraction
using KeldyshContraction: position

@qfields ϕ::Destroy(Classical) 

position(ϕ)
```

### Field Algebra

```@example API
using Term, KeldyshContraction # hide
Term.typestree(KeldyshContraction.QTerm) # hide
```

```@docs
KeldyshContraction.QTerm
KeldyshContraction.QMul
KeldyshContraction.QAdd
arguments(::KeldyshContraction.QMul)
arguments(::KeldyshContraction.QAdd)
```

The properties of the expression can be checked using:

```@docs
KeldyshContraction.is_bulk
KeldyshContraction.is_conserved
KeldyshContraction.is_physical
```

## Systems

```@docs
InteractionLagrangian
```

## Wick Contraction
  
```@docs
wick_contraction
```

### Propagator

```@docs
KeldyshContraction.PropagatorType
DressedPropagator
KeldyshContraction.matrix(::DressedPropagator)
```

### Self-Energy

```@docs
KeldyshContraction.SelfEnergy
KeldyshContraction.matrix(::KeldyshContraction.SelfEnergy)
```
