# Spectral/dispersive collision basis

The exact off-shell Kadanoff--Baym collision expression still contains retarded and advanced propagator factors. Before any quasiparticle or mass-shell approximation, these causal factors can be rewritten exactly in the spectral/dispersive basis

```math
D = \frac{G^R+G^A}{2},
\qquad
A = i(G^R-G^A),
```

so that

```math
G^R = D-\frac{i}{2}A,
\qquad
G^A = D+\frac{i}{2}A.
```

[`spectral_dispersive_collision`](@ref) applies these identities term by term to both affine pieces of an [`OffShellCollisionExpression`](@ref). Keldysh-derived spectral factors remain spectral and keep their statistical weight.

This transformation is exact. In particular, it does **not** replace `A` by a quasiparticle delta function and it does **not** replace `D` by a principal-value energy denominator. Those operations require an explicit dispersion relation and frequency elimination and belong to the subsequent projection stage.

The representation deliberately reuses the exact `KineticTerm` carrier for routed momenta, field families, equal-time regularisation shifts, topology and derivative-generated momentum polynomials. Original retarded/advanced line labels are normalized out of the carrier, while a fixed-size [`SpectralDispersiveKind`](@ref) vector records whether each physical factor is spectral or dispersive. This allows algebraically identical contributions coming from different R/A spellings to merge and cancel exactly.

A basic consistency identity is therefore represented literally:

```math
i(G^R-G^A)=A.
```

The two dispersive pieces cancel under duplicate merging, leaving one spectral factor with unit coefficient.

## Why this stage precedes shell projection

The separation matters for mixed elastic--inelastic sectors. A product involving a spectral factor and a dispersive factor can vanish only after the relevant frequency integration is shown to implement the appropriate Kramers--Kronig relation. Conversely, a surviving dispersive factor produces principal-value support rather than an energy-conserving delta function. Keeping `A` and `D` explicit prevents either behavior from being imposed by graph topology or by a heuristic simplifier.

Equal-time regularisation metadata also remains present at this stage. It must survive until the dedicated frequency/equal-time projection handles the regularised loss tadpole; no generic causal decomposition is allowed to erase it.

```@docs
SpectralDispersiveKind
SpectralDispersiveTerm
SpectralDispersiveExpression
SpectralDispersiveCollision
spectral_dispersive_collision
spectral_dispersive_kinds
```
