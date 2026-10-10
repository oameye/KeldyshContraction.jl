# Fourier and Wigner representations

`fourier_transform` converts coordinate-space diagrams to an exact momentum-space representation. Momentum conservation is solved symbolically and coordinate derivatives become momentum polynomials.

For a perturbative self-energy, Fourier routing must happen before amputation:

```julia
GF = fourier_transform(G)
ΣF = SelfEnergy(GF)
```

A coordinate-space `SelfEnergy` has already removed its external propagators and therefore cannot prove that derivative momentum factors at those endpoints were retained. Direct `fourier_transform(::SelfEnergy)` is rejected for that reason.

The 2PI route uses the same ordering. `twopi_fourier_self_energy` cuts the target full-propagator line while its interaction endpoints still exist, attaches the physical external propagator segments, routes and derivative-lowers the complete two-point graph, and only then amputates those segments.

```julia
ΣF = KeldyshContraction.twopi_fourier_self_energy(Γ2, ψ)
```

```@docs
fourier_transform
KeldyshContraction.twopi_fourier_self_energy
```

`wigner_transform` then separates centre and relative coordinates. The gradient order is explicit; the current kinetic compiler uses `Val(0)`.

```julia
ΣW = wigner_transform(ΣF; gradient_order=Val(0))
```

```@docs
wigner_transform
```

Both transformations preserve the semantic retarded, advanced, and Keldysh components and the physical target-field provenance. Finite equal-time regularisation is also preserved until the kinetic reduction resolves it.

Users normally do not need to construct or dispatch on the Fourier/Wigner intermediate types. Pass the returned object directly to the next stage of the workflow.

See [Quantum kinetic theory](../theory/kinetics.md) for the approximation made after the Wigner representation is formed.