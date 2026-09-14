# Fields and algebra

Declare physical field families with `@qfields`. A family is either bosonic or fermionic; its Keldysh component is selected by indexing.

```julia
using KeldyshContraction

@qfields ϕ::Boson ψ::Fermion
c, q = ϕ[Classical], ϕ[Quantum]
ψ1, ψ2 = ψ[One], ψ[Two]
```

```@docs
KeldyshContraction.@qfields
Boson
Fermion
FieldFamily
Field
field_family
```

`bar` changes field orientation. For fermions it constructs the independent barred Grassmann variable; it is not an operator adjoint.

```@docs
bar
```

Coordinate derivatives remain symbolic until Fourier transformation.

```julia
∂xψ = partial(ψ1, :x)
∂xyψ = partial(partial(ψ1, :x), :y)
```

```@docs
partial
derivatives
```

Bosonic products commute. Fermionic products obey Grassmann antisymmetry and nilpotency. Users normally work through ordinary `+`, `-`, `*`, and powers rather than constructing the package's algebra containers directly.

For finite-Trotter vertices, regularisation labels are attached by calling a field with `Regularisation.Plus` or `Regularisation.Minus`.

```@docs
KeldyshContraction.Regularisation
```

See [Equal-time regularisation](../theory/regularisation.md) for the physical meaning of these shifts.