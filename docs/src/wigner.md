# Wigner representation

The Wigner layer is a semantic representation built on the exact Fourier IR. It does not perform graph momentum routing, spectral reduction, distribution substitution, quasiparticle projection, or collision construction.

For a two-point object with coordinates ``x_1`` and ``x_2``, define center and relative coordinates

```math
X = \frac{x_1+x_2}{2},
\qquad
s = x_1-x_2,
```

and Fourier transform the relative coordinate ``s`` to its conjugate momentum ``p``. The Fourier convention and derivative-generated momentum factors are inherited unchanged from the Fourier layer.

## Initial homogeneous regime

The first production Wigner representation supports translation-invariant, homogeneous systems at gradient order zero:

```julia
ΣW = wigner_transform(Σk; gradient_order=Val(0))
```

At this order the numerical graph data are operationally identical to the exact Fourier representation, but the result is a distinct Wigner-space object with explicit semantics:

- the external relative-coordinate momentum is recorded as the Wigner momentum;
- exact loop-momentum routing is preserved;
- exact derivative/kinematic `MomentumPolynomial` factors are preserved without re-derivation;
- `HomogeneousWignerContext` records that no explicit center-coordinate dependence is represented;
- gradient order is visible statically in the Wigner result type.

The Wigner transformation is statistics-neutral. Bosonic and fermionic inputs use the same Wigner IR and transformation machinery; statistics-specific distribution identities belong downstream.

## Transform before coordinate-space amputation

Derivative-decorated interactions require the same provenance-preserving order established by the Fourier layer:

```julia
Gk = fourier_transform(G)
Σk = SelfEnergy(Gk)
ΣW = wigner_transform(Σk; gradient_order=Val(0))
```

A coordinate-space `SelfEnergy` may already have lost derivative metadata carried by amputated external propagators, so `wigner_transform(::SelfEnergy)` is rejected. A coordinate `DressedPropagator` may be passed directly as a convenience; it is Fourier-transformed first.

## Gradient-expansion seam

The public approximation parameter is explicit:

```julia
wigner_transform(x; gradient_order=Val(0))
```

Only `Val(0)` is currently supported. Any nonzero gradient order fails explicitly rather than silently returning the homogeneous result.

The IR is designed so future nonzero orders can represent the Moyal/star-product expansion

```math
A \star B
= A B
+ \frac{i}{2}
\left(
\partial_X A\,\partial_p B
- \partial_p A\,\partial_X B
\right)
+ \cdots.
```

No nonzero-gradient term is implemented in this layer. Likewise, spectral/statistical reduction and kinetic/collision approximations remain separate downstream stages.
