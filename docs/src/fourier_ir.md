# Fourier routing IR

The momentum/Fourier layer is deliberately separate from coordinate-space diagram identity.
`Diagram` continues to represent coordinate-space graph and field data. Exact momentum routing is
attached only when constructing a `FourierDiagram`.

For a two-point diagram, the routing convention follows the directed propagator orientation: momentum
flows from an edge's `in` endpoint to its `out` endpoint. The incoming external propagator therefore
carries the canonical external momentum from `In()` into the bulk, while the outgoing propagator carries
the same momentum from the bulk to `Out()`.

After removing those fixed external legs from the unknown internal routing problem, the exact system is

```math
B p = q_{\mathrm{ext}},
```

and is solved as an affine exact-rational system. The first momentum-basis variable is the external
momentum. Remaining basis variables are deterministic loop momenta selected by left-to-right exact RREF
free columns. No floating-point linear algebra or topology-specific sign heuristic is used.

`FourierDiagram` stores the original coordinate-space `Diagram`, the canonical momentum basis, and the
exact linear momentum assigned to every coordinate-space edge. Its equality and hash include the routing,
so physically different momentum-space objects cannot alias merely because their coordinate graphs are
equal.

## Derivative momentum factors

Coordinate derivatives are lowered only after exact routing. The Fourier convention is

```math
A(X,k)=\int ds\; e^{-i(\mathbf{k}\cdot\mathbf{s}-\epsilon t)}
A(X+s/2,X-s/2).
```

At homogeneous gradient order zero, the inverse transform therefore gives the endpoint rules

```math
\partial_{x_{\rm out}^j} \mapsto +i p_j,
\qquad
\partial_{x_{\rm in}^j} \mapsto -i p_j,
```

for spatial axes, while time derivatives carry the opposite sign because the Fourier phase is
`\mathbf{k}\cdot\mathbf{s}-\epsilon t`.

`MomentumComponent` stores one axis component of an exact `LinearMomentum`. Products are represented by
canonical `MomentumMonomial`s, and finite sums by `MomentumPolynomial{C}`. Fourier derivative lowering
uses exact `ComplexRationals` coefficients for powers of `i`; the ordinary numeric diagram coefficient
remains separate and is never widened to a symbolic momentum expression.

A routed `FourierDiagram{...,Nothing}` is converted by `lower_fourier_derivatives` to a
`FourierDiagram` whose `kinematic` field is a concrete exact momentum polynomial. The routed linear
momentum is retained as one factor rather than expanded into an external/loop expression tree. This is the
representation needed for derivative/p-wave vertices while preserving deterministic routing identity.

## Public Fourier transformation

`fourier_transform` is the public translation-invariant coordinate-to-momentum operation. It does not
perform a center/relative-coordinate Wigner transform or a gradient expansion.

For one `Diagram`, `fourier_transform` returns its canonically routed diagram with all coordinate
derivatives lowered to the exact kinematic polynomial. For `Diagrams`, the result is a
`FourierDiagrams` collection keyed by the derivative-consumed routed graph. Numerical diagram
coefficients and momentum-polynomial kinematics remain separate. Coordinate diagrams that differ only by
derivative placement can therefore meet under one physical Fourier graph while retaining distinct exact
kinematic contributions.

`fourier_transform(::DressedPropagator)` applies the same transformation component-wise and preserves the
coefficient representation, statistics, perturbation order, and static graph shape in
`FourierDressedPropagator`.

For derivative-decorated interactions, self-energy amputation must happen after Fourier lowering:

```julia
Gk = fourier_transform(G)
Σk = SelfEnergy(Gk)
```

This order is semantically required. A coordinate-space `SelfEnergy` has already removed its external
propagators, so it does not retain enough provenance to prove that no derivative metadata was discarded
at those endpoints. Direct `fourier_transform(::SelfEnergy)` is therefore rejected rather than silently
producing an incomplete momentum factor.

The Fourier-space self-energy retains the routed internal graph and the complete kinematic factor computed
before amputation. In particular, derivative factors carried by external legs remain present in the
self-energy contribution even though those propagators are no longer stored in the amputated graph.

The former experimental `wigner_transform` algorithm is not the Wigner layer of this architecture. The
actual center/relative-coordinate Wigner representation and gradient expansion are separate follow-up
work built on this Fourier IR.
