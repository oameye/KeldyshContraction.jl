# 2PI effective action

`TwoPIEffectiveAction` is the diagrammatic entry point for self-consistent two-particle-irreducible
(2PI) approximations.

At fixed interaction order it generates the connected vacuum skeleton contribution

```math
\Gamma_2[G]
```

directly from the interaction action. The intended compiler chain is

```text
interaction action
    -> TwoPIEffectiveAction
    -> graph functional derivative
    -> proper self-energy Σ[G]
    -> Dyson / Kadanoff--Baym representation
    -> Fourier / Wigner / kinetic compiler
```

Only the first stage is implemented here. `TwoPIEffectiveAction` does not solve a Dyson or
Kadanoff--Baym equation.

```@docs
TwoPIEffectiveAction
KeldyshContraction.twopi_terms
```

## Formal full propagators

The effective action is a functional of the **formal full propagator matrix**. This is different
from generating an ordinary physical vacuum expectation value.

In particular, KC's perturbative Wick path normally applies physical Keldysh identities while the
diagrams are generated: the quantum--quantum propagator is zero and causal vacuum loops vanish.
Those reductions must not be applied before differentiating ``\Gamma_2``. A contribution that
vanishes after restriction to the physical propagator manifold can still have a nonzero functional
derivative and therefore contribute to the proper self-energy.

For that reason the 2PI generator stores formal endpoint contractions directly and retains all
Keldysh matrix components, including `qq`, until a later differentiation/projection stage. It still
enforces field-family compatibility, barred/unbarred orientation, exact bosonic/fermionic Wick
parity, connectedness, and two-particle irreducibility.

## Fixed shape

The current constructor follows KC's explicit static-shape convention:

```julia
Γ2 = TwoPIEffectiveAction(L, Val(order), Val(edges))
```

For an interaction whose monomials each contain `N` fields,

```math
2E = O N,
```

where `O` is the interaction order and `E` is the number of full propagator lines in the vacuum
diagram.

The exact perturbative parameter monomial is retained. Semantic metadata are available through
`order`, `statistics`, `parameters`, `field_families`, and `twopi_terms`.

## Relation to `skeleton_self_energy`

`skeleton_self_energy` remains an independent perturbative oracle. It selects the 2PI subset of an
ordinary generated self-energy. It is useful for certifying the future graph-functional derivative,
but it is not the semantic foundation of the 2PI compiler.

The intended invariant is eventually

```text
functional_derivative(TwoPIEffectiveAction(...))
    == independently generated skeleton_self_energy(...)
```

diagram by diagram, with exact multiplicities and statistics signs.

## Current interaction boundary

This first slice accepts interactions already representable by `InteractionLagrangian`. Its
barred/unbarred conservation invariant is intentionally unchanged.

The cubic charged Hubbard--Stratonovich interaction used for two-body loss,

```math
\bar\chi\psi^2 + \chi\bar\psi^2,
```

requires an auxiliary/charged interaction representation because the auxiliary field carries a
different U(1) charge. That extension belongs to the next 2PI layers rather than weakening the
existing interaction contract.
