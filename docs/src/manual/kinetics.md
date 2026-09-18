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
LoopMomentumQuotientWorkspace
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

For the underlying equations and approximation boundaries, see [Quantum kinetic theory](../theory/kinetics.md).
