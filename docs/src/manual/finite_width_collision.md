# Finite-width collision reduction

The finite-width branch acts on the complete spectral/dispersive collision before the sharp-quasiparticle frequency limit. It preserves the exact affine Kadanoff--Baym structure

```math
I_{\mathrm{coll}}=I_0+F_{\mathrm{target}}(k)I_1
```

while classifying each routed term with the analytic Lorentzian backends. Independent or repeated spectral lines use the factorized reduction, the regular `L+1`-line / `L`-loop sector uses the residual Cauchy convolution, and structures outside those analytic domains remain explicitly unsupported.

This stage does not supply a linewidth model, choose the microscopic external spectral frequency, lower `F` to occupations, or integrate momenta. In particular, the current statistical factors remain the quasiparticle/statistical factors carried by the kinetic IR; this is not a general solver for a frequency-dependent distribution `F(ω,k)`.

```@docs
KeldyshContraction.FiniteWidthFrequencyReductionKind
KeldyshContraction.FiniteWidthFrequencyTerm
KeldyshContraction.FiniteWidthFrequencyExpression
KeldyshContraction.FiniteWidthFrequencyCollision
KeldyshContraction.finite_width_frequency_collision
```

## Spectral-model evaluation

After structural classification, an explicit `LorentzianSpectralModel` can evaluate every analytically resolved frequency factor. Factorized terms use the exact integrated Lorentzian powers, while residual-shell terms use the Cauchy convolution at an explicitly supplied microscopic spectral frequency `ω_external`.

Resolved and unsupported terms remain separate and both retain their original `FiniteWidthFrequencyTerm` key. Missing spectral data are therefore an error rather than an implicit sharp-shell limit or cutoff. This operation still does not lower statistical factors to occupations, integrate loop momenta, or introduce collective-mode dynamics; `ω_external` is not the center-time collective response frequency.

```@docs
KeldyshContraction.FiniteWidthEvaluatedExpression
KeldyshContraction.FiniteWidthEvaluatedCollision
KeldyshContraction.evaluate_finite_width_collision
```

## Occupation lowering at fixed spectral data

Once the spectral weights have been evaluated, the resolved branch can be lowered with the same statistics-generic algebra as the strict-quasiparticle collision compiler,

```math
F_S=1+2\sigma_S n_S,
\qquad
C_n=\frac{\sigma_S}{2}I_{\mathrm{package}}.
```

The original affine split is recombined here: only the distribution-coefficient branch receives the external target factor before `F -> n`. Unsupported finite-width frequency structures remain explicit and are not assigned an occupation-space value.

This stage treats the supplied `LorentzianSpectralModel` as fixed. If a linewidth is itself generated from the state, `Γ=Γ[n]`, a self-consistent linear response must also differentiate that spectral input. The resulting `δΓ` chain-rule contribution is a later response layer rather than part of fixed-model occupation lowering.

```@docs
KeldyshContraction.FiniteWidthOccupationCollision
KeldyshContraction.finite_width_occupation_collision
```
