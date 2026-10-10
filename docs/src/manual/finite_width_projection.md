# Finite-width response projection

A fully resolved finite-width collision response can be projected with the same explicit left/right basis and projection functional used by the strict-quasiparticle collision stack. The projection layer does not distinguish the microscopic origin of the response: both a fixed-`ω_external` `FiniteWidthCollisionLinearization` and an `ExternalSpectralCollisionLinearization`, in which the external spectral line has already been integrated, use the same matrix assembly.

For

```math
\delta C=\sum_{s,a} c_{s,a}\,\delta n_a,
```

a right-basis closure supplies the perturbation amplitude `\psi_j(a)` and the left functional performs the explicit downstream phase-space or moment projection,

```math
K_{ij}=\sum_{s,a}\mathcal P_i[s,a,c_{s,a}\psi_j(a)].
```

The coefficients `c_{s,a}` already contain the response appropriate to their microscopic boundary. In particular, an externally projected self-consistent response contains the complete product rule

```math
\delta(WP)=\bar W\,\delta P+\bar P\,\delta W
```

with the external spectral line included in `\delta W`. Projection therefore does not re-evaluate spectral data or introduce another linewidth approximation.

Projection is accepted only when the finite-width response has no unresolved offset or distribution sectors. The core still does not choose a trap model, phase-space quadrature, Chapman--Enskog closure, collective response frequency, or pole equation. At the current gradient order the projected operator is static; a frequency-dependent `M(\Omega)` requires a separate kinetic producer with genuine center-time memory.

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
