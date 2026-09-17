# Kinetic reduction

For a selected perturbative dressed-propagator sector `G`, the ordinary user API is one operation:

```julia
C = collision_kernel(G)
```

```@docs
collision_kernel
```

`collision_kernel(G)` runs the homogeneous zeroth-gradient kinetic compiler and returns the final strict-quasiparticle `CollisionKernel`. It is intentionally strict: genuine causal-frequency blockers and unresolved finite-Trotter states are reported rather than assigned an arbitrary finite value.

## Inspecting the derivation

The high-level call is exactly the composition below. These transformations remain public because they are useful for derivations, diagnostics, and research on intermediate approximations:

```julia
GF = fourier_transform(G)
ΣF = SelfEnergy(GF)
ΣW = wigner_transform(ΣF; gradient_order=Val(0))
KΣ = kinetic_expression(ΣW)
I = off_shell_collision_expression(KΣ)
SD = spectral_dispersive_collision(I)
R = reduce_frequency_collision(SD)
N = occupation_reduced_expression(R)
Q = quotient_loop_momenta(N)
C = collision_kernel(Q)
```

The canonical examples execute this chain one stage at a time. Each intermediate object has compact `text/plain` and `text/latex` displays so Documenter shows the represented physics rather than the implementation fields.

## Off-shell collision functional

`OffShellCollisionExpression` is the stable physical boundary before the strict-quasiparticle approximation. It stores the exact homogeneous Kadanoff--Baym collision functional affinely in the external statistical distribution,

```math
I_{\mathrm{coll}} = I_0 + F_{\mathrm{target}}(k) I_1,
\qquad I_0=i\Sigma^K,
\qquad I_1=-A_\Sigma.
```

No shell delta function, occupation substitution, principal-value reduction, or finite-width prescription has been applied at this stage. The semantic accessors expose the two affine pieces and the retained Wigner context without exposing internal term dictionaries.

```@docs
KeldyshContraction.OffShellCollisionExpression
KeldyshContraction.collision_offset
KeldyshContraction.collision_distribution_coefficient
KeldyshContraction.gradient_order
KeldyshContraction.wigner_context
```

## Finite-width spectral data

Finite-width analysis is a separate branch from the strict-quasiparticle reduction. The first analytic model is the normalized Lorentzian/Breit--Wigner line

```math
A(\omega)=\frac{Z\Gamma}{(\omega-E)^2+\Gamma^2/4}.
```

The line center, linewidth, and residue are explicit data. No global broadening convention is hidden in the compiler. For an on-shell retarded self-energy the package uses the explicit convention

```math
E_k=\varepsilon_k+\operatorname{Re}\Sigma^R_k,
\qquad
\Gamma_k=-2\operatorname{Im}\Sigma^R_k,
```

while the quasiparticle residue `Z` remains an independently supplied datum. Frequency-dependent self-consistency and derivative corrections to `Z` are deliberately outside this algebraic conversion.

For integer multiplicity `m`, the analytic line integral is

```math
\int\frac{d\omega}{2\pi}A(\omega)^m
=
\binom{2m-2}{m-1}\frac{Z^m}{\Gamma^{m-1}}.
```

In particular,

```math
\int\frac{d\omega}{2\pi}A(\omega)=Z,
\qquad
\int\frac{d\omega}{2\pi}A(\omega)^2=\frac{2Z^2}{\Gamma}.
```

`lorentzian_spectral_reduction` applies this identity only when an all-spectral term factorizes into independent full-rank line-frequency coordinates. Repeated lines are identified by exact physical `(field family, routed momentum)` identity. The operation always returns one concrete `LorentzianSpectralReduction`: resolved terms are tagged `LorentzianResolved`, while dispersive factors, finite Trotter shifts, rank-deficient routing, and nonfactorized structures are tagged `LorentzianUnsupported`. The analytic backend does not assign unsupported structures a cutoff value.

The regular `L + 1`-line / `L`-loop case is handled separately by `lorentzian_convolution_reduction`. The same exact frequency-routing matrix used by the strict quasiparticle reducer chooses a full-rank set of pivot lines and expresses the remaining line as

```math
\omega_d=\beta\,\omega_{\mathrm{ext}}+\sum_j c_j\omega_j.
```

Cauchy stability then gives an exact residual Lorentzian with

```math
\Gamma_{\mathrm{eff}}=\Gamma_d+\sum_j |c_j|\Gamma_j.
```

The stored residual `EnergyForm` is the same mismatch that becomes the strict energy shell. Consequently the zero-width limit recovers the strict `EnergyShell` with the same exact affine Jacobian and shell-normalization factor; no numerical limiting procedure is needed to certify this correspondence.

The resulting linewidth powers for repeated lines are stored separately from the perturbative `ParameterMonomial`. A nominal `γ²` term containing one repeated line therefore remains `γ²` with an explicit `Γ(q)^(-1)` factor until a linewidth scaling law is supplied.

```@docs
KeldyshContraction.SpectralLineIdentity
KeldyshContraction.LorentzianSpectralData
KeldyshContraction.LorentzianSpectralModel
KeldyshContraction.LorentzianSpectralFactor
KeldyshContraction.LorentzianSpectralWeight
KeldyshContraction.LorentzianSpectralReductionKind
KeldyshContraction.LorentzianSpectralReduction
KeldyshContraction.LorentzianConvolutionReductionKind
KeldyshContraction.LorentzianConvolutionReduction
KeldyshContraction.LinewidthPower
KeldyshContraction.WidthAwarePowerCounting
KeldyshContraction.spectral_line_family
KeldyshContraction.spectral_line_momentum
KeldyshContraction.spectral_data
KeldyshContraction.spectral_energy
KeldyshContraction.spectral_linewidth
KeldyshContraction.spectral_residue
KeldyshContraction.spectral_normalization
KeldyshContraction.spectral_squared_weight
KeldyshContraction.spectral_jacobian
KeldyshContraction.spectral_factors
KeldyshContraction.spectral_reduction_kind
KeldyshContraction.spectral_reduction_resolved
KeldyshContraction.convolution_reduction_kind
KeldyshContraction.convolution_reduction_resolved
KeldyshContraction.convolution_jacobian
KeldyshContraction.convolution_pivot_lines
KeldyshContraction.convolution_dependent_line
KeldyshContraction.convolution_coefficients
KeldyshContraction.convolution_external_coefficient
KeldyshContraction.convolution_energy_mismatch
KeldyshContraction.convolution_effective_linewidth
KeldyshContraction.spectral_line
KeldyshContraction.spectral_multiplicity
KeldyshContraction.linewidth_exponent
KeldyshContraction.linewidth_powers
KeldyshContraction.width_aware_power_counting
KeldyshContraction.lorentzian_spectral_data_from_retarded_self_energy
KeldyshContraction.lorentzian_spectral_value
KeldyshContraction.lorentzian_integrated_power
KeldyshContraction.lorentzian_spectral_reduction
KeldyshContraction.lorentzian_convolution_reduction
KeldyshContraction.evaluate_spectral_weight
KeldyshContraction.evaluate_spectral_convolution
```

## Spectral/statistical lowering

```@docs
kinetic_expression
off_shell_collision_expression
spectral_dispersive_collision
```

`off_shell_collision_expression` forms the complete Kadanoff--Baym collision identity before an on-shell approximation. `spectral_dispersive_collision` then separates spectral and causal/dispersive factors without discarding principal-value terms.

## Frequency and occupation reduction

```@docs
reduce_frequency_collision
occupation_reduced_expression
quotient_loop_momenta
```

The reduced frequency result separates finite shell/principal-value terms from unresolved equal-time or causal structures:

```@docs
KeldyshContraction.reduced_regular_terms
KeldyshContraction.reduced_blocked_terms
KeldyshContraction.reduced_trotter_terms
```

A blocked causal term is not silently assigned a finite strict-quasiparticle value. If a research calculation intentionally wants only the finite regular branch, that choice is explicit:

```julia
N = occupation_reduced_expression(R)
C_regular = collision_kernel(N)
```

## Final collision kernel

```@docs
CollisionKernel
KeldyshContraction.collision_kernel_terms
```

The kernel stores the occupation polynomial together with its exact kinematic and frequency support. For programmatic inspection use semantic accessors rather than internal representation types:

```@docs
KeldyshContraction.kinematic_factor
KeldyshContraction.frequency_support
```

## Moment and linear-response projection

The collision kernel can be linearized exactly before any phase-space closure. If

```math
C[n]=\sum_m c_m\prod_a n_a^{p_{ma}},
```

then `occupation_linearization` constructs the canonical first Fréchet derivative

```math
\delta C
=
\sum_a \delta n_a\,\frac{\partial C}{\partial n_a},
```

with exact integer multiplicities and exact polynomial coefficients. `linearize_collision_kernel` applies this operation sector by sector while leaving perturbative provenance, routed momenta, derivative kinematics, shell/PV support, gradient order and Wigner context unchanged.

A supplied background is a separate operation. `OccupationBackground(V, f)` declares both the evaluator and its value type. Evaluating the symbolic derivative at `\bar n` gives

```math
\delta C\big|_{\bar n}
=
\sum_a \delta n_a
\left.\frac{\partial C}{\partial n_a}\right|_{\bar n},
```

while retaining the exact variation atom and collision sector. Vanishing channels and sectors are removed canonically. The background evaluator may return numerical or symbolic numbers; it does not choose a phase-space measure, closure, or trap model.

Moment projection remains separate from phase-space integration. `CollisionMomentProjection` attaches an explicit test-function descriptor to nonlinear, symbolic-linearized, or background-evaluated collision data. The package provides exact built-in descriptors for

```math
\dot N = \int C,
\qquad
\dot E = \int \varepsilon_k C,
```

through `NumberMoment` and `EnergyMoment`. Arbitrary concrete user descriptors can be attached with `project_collision_moment`; downstream trap, measure, closure or quadrature code may interpret them without forcing that physics into the collision compiler. In particular, the later breathing-mode basis can use the same projection boundary without hard-coding `r²`, `r\!\cdot\!k`, or `k²` into the core IR.

```@docs
KeldyshContraction.OccupationLinearization
KeldyshContraction.LinearizedCollisionKernel
KeldyshContraction.CollisionMomentProjection
KeldyshContraction.NumberMoment
KeldyshContraction.EnergyMoment
KeldyshContraction.OccupationBackground
KeldyshContraction.BackgroundOccupationLinearization
KeldyshContraction.BackgroundLinearizedCollisionKernel
KeldyshContraction.occupation_linearization
KeldyshContraction.occupation_linearization_terms
KeldyshContraction.linearize_collision_kernel
KeldyshContraction.linearized_collision_terms
KeldyshContraction.occupation_background_value
KeldyshContraction.evaluate_occupation_polynomial
KeldyshContraction.evaluate_occupation_linearization
KeldyshContraction.background_occupation_linearization_terms
KeldyshContraction.evaluate_collision_background
KeldyshContraction.background_linearized_collision_terms
KeldyshContraction.project_collision_moment
KeldyshContraction.number_moment_projection
KeldyshContraction.energy_moment_projection
KeldyshContraction.moment_test_function
KeldyshContraction.projected_collision
KeldyshContraction.moment_weight
KeldyshContraction.linearize_collision_moment
```

This layer now supplies the exact linearized response at a user-provided background. It still does not evaluate phase-space integrals, close a moment basis, construct the projected matrix, or solve the collective-mode pole equation; those are the next #342 stages.

For the underlying equations and approximation boundaries, see [Quantum kinetic theory](../theory/kinetics.md).
