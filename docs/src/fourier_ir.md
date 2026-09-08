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

This layer intentionally does not yet lower coordinate derivatives to momentum polynomials. The
`kinematic` slot is therefore `nothing`; the following derivative/Fourier layer will populate it without
widening ordinary numeric diagram coefficients.
