# Exceptional frequency classification

The generic quasiparticle frequency reducer handles terms whose spectral constraints span all loop frequencies and whose residual shell/principal-value supports remain nonsingular. Some physically important terms sit outside that domain. They must be classified structurally rather than forced through the generic rule.

The internal exceptional-frequency classifier therefore returns a concrete classification record carrying the outcome kind, the original routed term, its spectral rank and loop count, and explicit proof provenance where available. This classifier is implementation machinery rather than part of the stable public API; users should consume the public reduced-frequency result types and accessors instead.

The current internal outcomes are:

- `FrequencyKramersKronigZero`: a conservative causal-orthogonality proof exists;
- `FrequencyTrotterRequired`: at least one equal-time line retains a nonzero regularisation shift and must be handled by the dedicated Trotter rule;
- `FrequencyUnresolved`: no implemented exceptional identity is strong enough to reduce the term.

Rank deficiency alone has no physical interpretation. In particular it does **not** mean that a term vanishes and it does **not** authorize replacing the remaining frequency structure by a principal-value denominator.

## Conservative Kramers--Kronig proof

A Kramers--Kronig zero is accepted only when one loop-frequency variable appears in exactly two propagator factors and those two factors are one spectral factor and one dispersive factor with the same causal denominator. All remaining propagators are then independent of that integration variable. Within the quasiparticle statistical convention the distribution weight on a Keldysh-derived spectral line depends on momentum, not on the integrated frequency, and therefore factors out as well. The remaining frequency integral is the causal orthogonality relation

```math
\int d\omega\, A(\omega) D(\omega) = 0,
\qquad
D = \frac{G^R+G^A}{2},
\quad
A = i(G^R-G^A).
```

For this outcome, the internal classification stores the spectral/dispersive witness lines and the isolated loop-frequency basis index. The proof deliberately requires the same family, endpoints, regularisation shift and exact routed momentum for the two causal factors. More general Hilbert-transform identities can be added later only with equally explicit provenance.

This rule is the structural mechanism needed for the analytically known vanishing mixed `gγ` three-multiplicity sector. It is intentionally narrower than that final physical acceptance test: if the generated `gγ` term does not expose this exact witness yet, it remains unresolved until the required algebraic identity is represented explicitly rather than being matched by topology or perturbative order.

## Equal-time terms

Any nonzero regularisation shift takes precedence over the ordinary causal classification. Such a term is internally classified as `FrequencyTrotterRequired`; no Plemelj--Sokhotski or Kramers--Kronig limit is applied. This preserves the factor-of-two information required by the regularised two-body-loss tadpole sector.
