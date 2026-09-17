# Moment and linear-response projection

The projection layer sits downstream of the canonical collision compiler. It does not change the
collision IR and does not assume a phase-space measure, trap model, quadrature rule, or moment
closure.

After exact occupation-space linearization and evaluation at a supplied background,

```math
\delta C = \sum_{s,a} c_{s,a}\,\delta n_a,
```

`CollisionProjectionBasis` stores two explicit bases:

- a left basis of observables or test functions ``\phi_i``;
- a right basis of perturbations ``\psi_j``.

The bases are deliberately distinct. The core does not assume a self-dual moment basis.

`CollisionPerturbationClosure` supplies the right-basis amplitudes ``\psi_j(a)`` on the exact
occupation variation channels. `CollisionProjectionFunctional` supplies the left projection of one
closed response channel. The projected collision matrix is therefore

```math
K_{ij}
=\sum_{s,a}
\mathcal P_i\!\left[s,a,c_{s,a}\psi_j(a)\right].
```

This is the finite-dimensional collision contribution required by the later collective-mode
problem. Analytic phase-space integration, trapped-gas closure, and the streaming part of the
kinetic operator remain explicit downstream physics. The representation also does not require a
future finite-width projected operator to be frequency independent; a later backend may construct
``M(\Omega)`` using the same explicit left/right basis semantics.

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
