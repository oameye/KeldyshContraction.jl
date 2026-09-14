# Propagators and self-energy

`DressedPropagator` generates the perturbative two-point function at a chosen interaction order. The number of propagator edges is static and is therefore supplied as a `Val`.

```julia
G = DressedPropagator(L, Val(2), Val(5))
```

Use `preserve_regularisation=true` when the interaction contains finite-Trotter shifts.

```@docs
DressedPropagator
```

The retarded, advanced, and Keldysh components are accessed semantically rather than through storage fields.

```@docs
KeldyshContraction.retarded_component
KeldyshContraction.advanced_component
KeldyshContraction.keldysh_component
KeldyshContraction.matrix
KeldyshContraction.order
KeldyshContraction.statistics
```

For direct diagrammatic work, `wick_contraction` exposes the underlying Wick expansion and `topologies` groups diagrams by graph topology.

```@docs
wick_contraction
topologies
```

`SelfEnergy` amputates the external lines and retains the irreducible two-point diagrams,

```julia
Σ = SelfEnergy(G)
```

with the same retarded/advanced/Keldysh accessors.

```@docs
SelfEnergy
```

For a sum of interactions, `DressedPropagator` and `SelfEnergy` preserve the perturbative parameter sectors. Index the result with the same symbolic parameter expressions used to label the interactions.

The package constructs the fixed-order 1PI self-energy; a nonlinear Dyson iteration is not performed automatically. See [Perturbation theory and self-consistency](../theory/perturbation.md).