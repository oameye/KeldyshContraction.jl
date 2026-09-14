# Fourier and Wigner representations

`fourier_transform` converts coordinate-space diagrams to an exact momentum-space representation. Momentum conservation is solved symbolically and coordinate derivatives become momentum polynomials.

```julia
ΣF = fourier_transform(Σ)
```

```@docs
fourier_transform
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