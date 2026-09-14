# Kinetic reduction

Starting from a Wigner-space self-energy, the collision kernel is obtained by composing the public transformations below.

```julia
KΣ = kinetic_expression(ΣW)
I = off_shell_collision_expression(KΣ)
SD = spectral_dispersive_collision(I)
R = reduce_frequency_collision(SD)
N = occupation_reduced_expression(R)
Q = quotient_loop_momenta(N)
C = collision_kernel(Q)
```

Each call has one job; keeping the stages explicit makes the approximation boundary inspectable.

## Spectral/statistical lowering

```@docs
kinetic_expression
off_shell_collision_expression
spectral_dispersive_collision
```

`off_shell_collision_expression` forms the Kadanoff--Baym collision identity before an on-shell approximation. `spectral_dispersive_collision` separates spectral and causal/dispersive factors without discarding principal-value terms.

## Frequency and occupation reduction

```@docs
reduce_frequency_collision
occupation_reduced_expression
quotient_loop_momenta
```

The reduced frequency result separates finite shell/principal-value terms from unresolved equal-time or causal structures. The following inspection functions are stable when that distinction matters:

```@docs
KeldyshContraction.reduced_regular_terms
KeldyshContraction.reduced_blocked_terms
KeldyshContraction.reduced_trotter_terms
```

A blocked causal term is not silently assigned a finite strict-quasiparticle value. In particular, genuine pinch singularities remain explicit.

## Final collision kernel

```@docs
CollisionKernel
collision_kernel
KeldyshContraction.collision_kernel_terms
```

The kernel stores the occupation polynomial together with its kinematic and frequency support. For programmatic inspection use the semantic accessors rather than internal storage:

```@docs
KeldyshContraction.kinematic_factor
KeldyshContraction.frequency_support
KeldyshContraction.parameters
```

For the underlying equations and approximations, see [Quantum kinetic theory](../theory/kinetics.md).