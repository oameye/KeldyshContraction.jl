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

Coordinate derivatives are lowered only after exact routing. The Fourier/Wigner convention is

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
