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
KeldyshContraction.Position
```

#### Field Constructors

It is expected to create the fields using the `Create` and `Destroy` constructors, together with the `@qfields` macro.

```@docs
KeldyshContraction.@qfields
```

The created fields are callable to change a property of the individual fields:

```@example API
using KeldyshContraction
using KeldyshContraction: In, Classical, position

@qfields ϕ::Destroy(Classical) 

position(ϕ(In))
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
```

The properties of the expression can be checked using:

```@docs
KeldyshContraction.isbulk
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
KeldyshContraction.Propagator
KeldyshContraction.PropagatorType
KeldyshContraction.propagator
DressedPropagator
```

### Self-Energy

```@docs
KeldyshContraction.SelfEnergy
```
