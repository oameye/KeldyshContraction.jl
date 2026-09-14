# Perturbation theory and self-consistency

Split the contour action into a Gaussian part and an interaction,

```math
S=S_0+S_{\mathrm{int}}.
```

Expanding `\exp(iS_{\mathrm{int}})` and applying Wick's theorem expresses interacting correlators in products of Gaussian two-point functions. For the two-point function this produces the usual diagrammatic series.

## From Wick contractions to the self-energy

The perturbative Green function may be reorganised into one-particle-irreducible two-point insertions,

```math
G=G_0+G_0\circ\Sigma\circ G_0
  +G_0\circ\Sigma\circ G_0\circ\Sigma\circ G_0+\cdots,
```

or equivalently

```math
G^{-1}=G_0^{-1}-\Sigma.
```

The self-energy is the 1PI two-point kernel. At fixed order,

```math
\Sigma=\Sigma^{(1)}+\Sigma^{(2)}+\cdots.
```

KeldyshContraction.jl performs the Wick expansion, keeps the exact combinatorial prefactors and fermionic permutation signs, removes reducible two-point diagrams, and returns the retarded, advanced, and Keldysh components of `\Sigma`.

Coordinate derivatives remain attached to their contraction endpoints at this stage. They become momentum polynomials only after Fourier transformation.

## Fixed-order and self-consistent approximations

Two approximations that are often conflated should be kept distinct.

A fixed-order calculation evaluates the retained diagrams with the Gaussian propagator,

```math
\Sigma\approx\Sigma^{(m)}[G_0].
```

A self-consistent or skeleton approximation keeps a chosen diagram class but dresses its internal lines,

```math
\Sigma\approx\Sigma^{(m)}[G],
\qquad
G^{-1}=G_0^{-1}-\Sigma[G].
```

The first statement selects perturbative diagrams; the second adds a nonlinear Dyson closure and therefore resums infinitely many ordinary perturbative contributions.

KeldyshContraction.jl currently constructs the explicit fixed-order diagrammatic self-energy. It does **not** iterate the Dyson equation to self-consistency. The generated 1PI expressions can be used as the diagrammatic input to a separate self-consistent solver, but that numerical closure is outside the present package.

## Multiple processes

Different interactions may be combined while keeping their perturbative parameters distinct. If

```math
S_{\mathrm{int}}=gS_g+\gamma S_\gamma,
```

then second order contains separate `g^2`, `g\gamma`, and `\gamma^2` sectors. The package preserves that provenance through the self-energy and kinetic pipeline rather than expanding all processes into one anonymous expression.

This becomes important for open systems: coherent and dissipative vertices can contribute different shell, principal-value, and singular sectors even when they share the same fields.

See [Equal-time regularisation](regularisation.md) for the additional ordering information required by coincident-time dissipative vertices.