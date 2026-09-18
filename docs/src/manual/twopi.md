# 2PI effective-action compiler

KeldyshContraction can organize self-consistent resummations at the level of the two-particle
irreducible (2PI) effective action rather than by filtering perturbative self-energy diagrams after
generation.

The intended compiler chain is

```text
interaction action
    -> connected 2PI vacuum skeletons Γ₂[G]
    -> exact graph functional derivative δΓ₂/δG
    -> proper self-energy Σ[G]
    -> Dyson/Kadanoff--Baym representation
    -> Fourier/Wigner/kinetic compiler
```

The first implementation slice is deliberately diagrammatic. It constructs exact connected vacuum
skeletons from the existing Wick algebra and preserves exact coefficients, field-family provenance,
statistics, interaction order, and parameter monomials. It does not solve self-consistent Dyson or
Kadanoff--Baym equations.

The existing `skeleton_self_energy` selector remains useful as an independent perturbative oracle:
a self-energy obtained by differentiating `Γ₂` should agree diagram by diagram with the corresponding
2PI subset of the ordinary perturbative self-energy. It is not the semantic foundation of the 2PI
compiler.

A later interaction layer will support charged auxiliary fields such as the cubic Hubbard--Stratonovich
vertex `χ̄ ψ² + χ ψ̄²`. The current `InteractionLagrangian` intentionally retains its existing
barred/unbarred conservation invariant and is not weakened for that purpose.
