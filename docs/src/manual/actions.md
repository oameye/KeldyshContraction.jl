# Actions and interactions

`InteractionLagrangian` turns a symbolic field expression into a perturbative interaction vertex. The expression must be conserved, physical, and contain both Keldysh components.

```julia
@qfields ϕ::Boson
c, q = ϕ[Classical], ϕ[Quantum]

elastic = -(1 // 2) * ((c^2 + q^2) * bar(c) * bar(q) + c * q * (bar(c)^2 + bar(q)^2))
L = InteractionLagrangian(elastic, :g)
```

```@docs
InteractionLagrangian
```

The second argument labels the perturbative parameter. Distinct physical processes can be added directly,

```julia
L = InteractionLagrangian(elastic, :g) + InteractionLagrangian(loss, :γ)
```

and remain separated into `g`, `γ`, `g^2`, `gγ`, and `γ^2` sectors at higher order.

For interactions involving several physical field families, select the external family with the `target` keyword when constructing a propagator. The stable provenance accessors are:

```@docs
KeldyshContraction.field_families
KeldyshContraction.target_family
KeldyshContraction.parameters
```

The ordinary user should not manipulate `LagrangianSum` or perturbative bookkeeping containers directly; addition, indexing by parameter expressions, and the high-level propagator constructors provide the supported interface.

See [Perturbation theory and self-consistency](../theory/perturbation.md) for the diagrammatic meaning of the perturbative sectors.