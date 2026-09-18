# 2PI effective action

`TwoPIEffectiveAction` is the diagrammatic entry point for self-consistent two-particle-irreducible
(2PI) approximations.

At fixed interaction order it generates the connected vacuum skeleton contribution

```math
\Gamma_2[G]
```

directly from the interaction action. The current compiler milestone is

```text
interaction action
    -> TwoPIEffectiveAction
    -> graph functional derivative
    -> physical SelfEnergy
```

This reaches self-energy compiler parity with KC's existing 1PI route. Dyson/Kadanoff--Baym
equation IR, self-consistent Fourier/Wigner routing, and numerical solution are deliberately
deferred beyond this milestone.

```@docs
TwoPIEffectiveAction
ChargedInteractionLagrangian
KeldyshContraction.field_charge
KeldyshContraction.twopi_terms
KeldyshContraction.twopi_self_energy
KeldyshContraction.twopi_self_energy_terms
KeldyshContraction.physical_self_energy
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
Keldysh matrix components, including `qq`, until functional differentiation has been performed. It
still enforces field-family compatibility, barred/unbarred orientation, exact bosonic/fermionic
Wick parity, connectedness, and two-particle irreducibility.

## Functional derivative and physical projection

KC represents a complex bosonic propagator as one oriented line from an unbarred field to its
barred partner. In this non-Nambu representation the stationary 2PI functional uses an `i Tr`
propagator trace, and `twopi_self_energy` implements

```math
\Sigma_{a b}(x,x') = i\,\frac{\delta\Gamma_2}{\delta G_{b a}(x',x)}.
```

This is equivalent to the familiar

```math
\Sigma = 2i\,\frac{\delta\Gamma_2}{\delta G^{T}}
```

when a complex field is written in an explicitly doubled Nambu representation: the Nambu trace
then carries an additional factor `1/2`. The two-body-loss 2PI derivation uses that doubled
convention, while KC's graph IR stores each oriented complex propagator only once. Importing the
Nambu `2i` prefactor directly into the oriented graph representation would therefore double-count
the self-energy.

The derivative is an exact graph operation. Each target-family full-propagator line is cut once,
the two cut vertices remain distinguished during canonicalization, and equivalent cuts accumulate
with their exact multiplicity.

The result remains a formal self-energy functional. Internal `qq` propagators are not removed at
this stage. `physical_self_energy` performs the subsequent restriction to the physical bosonic
Keldysh manifold and maps the surviving channels to the existing retarded, advanced, and Keldysh
`SelfEnergy` representation.

This order matters: cutting a formal `G_qq` line is what generates the physical Keldysh
self-energy channel. Imposing `G_qq = 0` before differentiation would erase it.

The graph-cut machinery is statistics-generic, but the derivative normalization is convention
dependent. The present implementation certifies the oriented complex-boson factor `i` against the
independently generated skeleton self-energy. Fermionic 2PI differentiation remains explicitly
unsupported until its convention is independently derived and tested.

## Charged auxiliary interactions

`ChargedInteractionLagrangian` extends the interaction boundary without weakening the existing
`InteractionLagrangian` conservation contract. Each field family receives an immutable integer
U(1) charge, and every monomial must have zero orientation-weighted charge.

For the Hubbard--Stratonovich two-body-loss vertex,

```math
\bar\chi\psi^2 + \chi\bar\psi^2,
```

the assignments are `qψ = 1` and `qχ = 2`. Thus `bar(χ) ψ ψ` is neutral because
`-2 + 1 + 1 = 0`, although it does not contain equal raw numbers of barred and unbarred fields.
The ordinary interaction representation continues to reject this monomial.

The formal vacuum matcher evaluates the complete selected vertex product before pairing. Products
that do not contain equal total numbers of unbarred and barred fields have no complete Wick
matching and contribute zero. At leading nonzero order the charged cubic interaction therefore
selects one vertex of each orientation and produces the connected `G²D` 2PI skeleton without an
HS-specific topology path.

## Self-energy solver parity

The high-level 2PI route lowers to the same physical `SelfEnergy` compiler object as the existing
1PI route:

```julia
Σ1 = SelfEnergy(G)          # existing 1PI compiler
Σ2 = SelfEnergy(Γ2)         # single-family 2PI compiler
Σψ = SelfEnergy(Γ2, ψ)      # explicit target in a multi-family theory
Ωχ = SelfEnergy(Γ2, χ)
```

`SelfEnergy(Γ2, target)` is only a thin semantic lowering of

```julia
physical_self_energy(twopi_self_energy(Γ2, target))
```

and does not introduce another graph derivative or physical projection path. The target-free form
is accepted only when `Γ2` contains one field family; multi-family theories require an explicit
target.

Once lowered, `order`, `statistics`, `parameters`, `target_family`, `matrix`, and the retarded,
advanced, and Keldysh component accessors have exactly the same meaning as for the existing 1PI
`SelfEnergy` result. This `SelfEnergy` boundary is the stopping point of the current 2PI
infrastructure program.

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
ordinary generated self-energy; it is not the semantic foundation of the 2PI compiler.

The certification invariant is

```text
SelfEnergy(TwoPIEffectiveAction(...))
    == independently generated skeleton_self_energy(...)
```

diagram by diagram, with exact multiplicities and Keldysh components.
