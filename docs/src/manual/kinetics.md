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

### Generated linewidth kernel

The linewidth can also be generated from the same microscopic self-energy. `spectral_self_energy_kernel` starts from a `KineticSelfEnergy`, forms

```math
A_\Sigma=i(\Sigma^R-\Sigma^A),
```

runs that expression through the same exact Trotter/frequency reduction used by the collision compiler, and substitutes `F = 1 + 2σn`. The important normalization distinction is explicit: this path uses `occupation_substitute` directly and never applies the collision conversion `C_n=σ I/2`.

The resulting `SpectralSelfEnergyKernel` therefore represents the quasiparticle linewidth object itself. Finite shell/PV sectors are stored as occupation polynomials, while causal blockers and unresolved finite-Trotter states remain inspectable. No background distribution, momentum quadrature, trap model, or finite-width prescription is inserted by this step.

For bosonic two-body loss, the first-order generated oracle is

```math
A_{\Sigma,\gamma}(k)=\Gamma_\gamma(k)
=4\gamma\int_q n_q,
```

with the vacuum contribution cancelled by the equal-time regularisation. Evaluating this kernel on a physical background supplies the `Γ` entering the Lorentzian spectral data; that downstream background evaluation is separate from the compiler object itself.

```@docs
KeldyshContraction.SpectralSelfEnergySector
KeldyshContraction.SpectralSelfEnergyKernel
KeldyshContraction.spectral_self_energy_terms
KeldyshContraction.spectral_self_energy_blocked_terms
KeldyshContraction.spectral_self_energy_trotter_terms
KeldyshContraction.spectral_self_energy_kernel
```

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

### Supplied occupation backgrounds

A supplied occupation state is represented by a typed evaluator. This is deliberately only an algebraic substitution boundary: it evaluates exact occupation polynomials but does not choose a momentum integral, trap, quadrature, or closure.

```@docs
KeldyshContraction.OccupationBackground
KeldyshContraction.occupation_background_value
KeldyshContraction.evaluate_occupation_polynomial
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

For the underlying equations and approximation boundaries, see [Quantum kinetic theory](../theory/kinetics.md).
