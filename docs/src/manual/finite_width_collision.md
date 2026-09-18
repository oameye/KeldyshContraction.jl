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

This stage treats the supplied `LorentzianSpectralModel` as fixed. If a linewidth is itself generated from the state, `\Gamma=\Gamma[n]`, a self-consistent linear response must also differentiate that spectral input.

```@docs
KeldyshContraction.FiniteWidthOccupationCollision
KeldyshContraction.finite_width_occupation_collision
```

## Spectral-data response

The finite-width spectral model has an explicit line-resolved first variation

```math
(E,\Gamma,Z)\longrightarrow (\delta E,\delta\Gamma,\delta Z).
```

For one line with multiplicity `m`,

```math
B_m(Z,\Gamma)
=\binom{2m-2}{m-1}\frac{Z^m}{\Gamma^{m-1}},
```

and the implementation differentiates this expression directly. In particular,

```math
\delta\!\left(\frac{2Z^2}{\Gamma}\right)
=\frac{4Z}{\Gamma}\,\delta Z
-\frac{2Z^2}{\Gamma^2}\,\delta\Gamma.
```

For residual Cauchy convolutions, the same exact stored frequency geometry is reused while propagating the line-energy, linewidth, and residue variations through the mismatch, effective linewidth, and residue product. The microscopic `ω_external` is held fixed. It remains distinct from any collective center-time response frequency.

A model variation is explicit: when a response channel is supplied, missing required line variations are errors rather than implicit zeros.

```@docs
KeldyshContraction.LorentzianSpectralDataVariation
KeldyshContraction.LorentzianSpectralModelVariation
KeldyshContraction.spectral_data_variation
KeldyshContraction.spectral_energy_variation
KeldyshContraction.spectral_linewidth_variation
KeldyshContraction.spectral_residue_variation
KeldyshContraction.lorentzian_integrated_power_variation
KeldyshContraction.evaluate_spectral_weight_variation
KeldyshContraction.evaluate_spectral_convolution_variation
```

## Self-consistent finite-width response

When the spectral data depend on the kinetic state, the collision response must differentiate both the statistical polynomial and its finite-width spectral weight. For each resolved physical term,

```math
\delta(WP)=\bar W\,\delta P+\bar P\,\delta W.
```

`OccupationBackground` evaluates the occupation polynomial and its exact Fréchet derivative at a supplied background. `OccupationSpectralResponse` separately maps explicit occupation-variation channels to line-resolved spectral-model variations. `FiniteWidthCollisionLinearization` retains the two contributions independently as well as their sum, with the original `FiniteWidthFrequencyTerm` as the key.

This separation is important for a self-consistent linewidth such as `\Gamma[n]`: the second term contains the response of the spectral width and must not be lost by linearizing only the fixed-model occupation polynomial. Unsupported finite-width structures remain explicit. The operation does not integrate loop momenta or introduce trap, Chapman--Enskog, moment-projection, or collective-mode semantics.

```@docs
KeldyshContraction.OccupationBackground
KeldyshContraction.occupation_background_value
KeldyshContraction.evaluate_occupation_polynomial
KeldyshContraction.OccupationLinearization
KeldyshContraction.occupation_linearization
KeldyshContraction.occupation_linearization_terms
KeldyshContraction.BackgroundOccupationLinearization
KeldyshContraction.evaluate_occupation_linearization
KeldyshContraction.background_occupation_linearization_terms
KeldyshContraction.OccupationSpectralResponse
KeldyshContraction.spectral_occupation_response_terms
KeldyshContraction.FiniteWidthCollisionLinearization
KeldyshContraction.finite_width_linearized_terms
KeldyshContraction.finite_width_occupation_response_terms
KeldyshContraction.finite_width_spectral_response_terms
KeldyshContraction.linearize_finite_width_collision
```

## Finite-dimensional projection

A fully resolved `FiniteWidthCollisionLinearization` can be projected with the same response-agnostic basis, perturbation closure, and left projection functional used by the strict-quasiparticle collision stack. For

```math
\delta C=\sum_{s,a} c_{s,a}\,\delta n_a,
```

a right-basis closure supplies the perturbation amplitude `\psi_j(a)` and the left functional performs the explicit downstream phase-space or moment projection,

```math
K_{ij}=\sum_{s,a}\mathcal P_i[s,a,c_{s,a}\psi_j(a)].
```

The finite-width coefficients `c_{s,a}` already contain the complete self-consistent product-rule response `\bar W\,\delta P+\bar P\,\delta W`. KC still does not choose the trap basis, phase-space quadrature, Chapman--Enskog closure, collective-mode frequency, or pole condition.

Projection requires complete finite-width resolution. If either the offset or distribution branch still contains unsupported sectors, `projected_collision_matrix` throws instead of silently discarding them.

```@docs
KeldyshContraction.CollisionProjectionBasis
KeldyshContraction.CollisionPerturbationClosure
KeldyshContraction.CollisionProjectionFunctional
KeldyshContraction.ProjectedCollisionMatrix
KeldyshContraction.left_projection_basis
KeldyshContraction.right_perturbation_basis
KeldyshContraction.perturbation_amplitude
KeldyshContraction.project_collision_channel
KeldyshContraction.projected_collision_matrix
```

## Generated microscopic linewidth

The Lorentzian linewidth can be produced from the same generated Keldysh self-energy rather than inferred from the collision kernel. `spectral_self_energy_kernel` forms

```math
A_\Sigma=i(\Sigma^R-\Sigma^A)=-2\operatorname{Im}\Sigma^R
```

and reduces it with the same exact causal-frequency and occupation algebra as the collision compiler, but without the collision normalization factor. The resulting `SpectralSelfEnergyKernel` retains its routed loop momentum, perturbative parameter, kinematic support, and unresolved causal/Trotter provenance.

Evaluating that kernel on an occupation background still leaves the loop integral explicit. For first-order bosonic two-body loss,

```math
\Gamma_\gamma(k)=4\gamma\int_q n_q,
```

so a downstream calculation may supply both

```math
\bar\Gamma=4\gamma\int_q\bar n_q,
\qquad
\delta\Gamma=4\gamma\int_q\delta n_q
```

and use them as the linewidth and linewidth variation in the generic finite-width response above. KC deliberately does not turn this boundary into a quadrature, trap, or mode-projection API.

```@docs
KeldyshContraction.SpectralSelfEnergySector
KeldyshContraction.SpectralSelfEnergyKernel
KeldyshContraction.BackgroundSpectralSelfEnergyKernel
KeldyshContraction.spectral_self_energy_terms
KeldyshContraction.background_spectral_self_energy_terms
KeldyshContraction.spectral_self_energy_blocked_terms
KeldyshContraction.spectral_self_energy_trotter_terms
KeldyshContraction.spectral_self_energy_kernel
KeldyshContraction.evaluate_spectral_self_energy_background
```
