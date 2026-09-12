# Public API policy

KeldyshContraction.jl distinguishes three levels of package surface.

1. **Exported workflow API** is imported by `using KeldyshContraction` and is intended for the ordinary symbolic-to-kinetic user workflow.
2. **Qualified public API** is marked with `SciMLPublic.@public`. These names are stable and documented, but users access them as `KeldyshContraction.name` (or import them explicitly) so the package does not pollute the caller namespace with representation and inspection utilities.
3. **Internal implementation details** are neither exported nor marked public. Their names, storage, and dispatch structure may change without constituting an API change.

On Julia 1.11 and newer the qualified public declarations use Julia's native `public` semantics. SciMLPublic provides the compatibility layer required by the package's Julia 1.10 support.

The intended high-level workflow is

```text
Field / partial / bar
  -> InteractionLagrangian
  -> DressedPropagator
  -> SelfEnergy
  -> fourier_transform
  -> wigner_transform
  -> kinetic_expression
  -> off_shell_collision_expression
  -> spectral_dispersive_collision
  -> reduce_frequency_collision
  -> occupation_reduced_expression
  -> quotient_loop_momenta
  -> collision_kernel
```

The intermediate stages remain explicit because they correspond to distinct mathematical operations and approximation boundaries. Public inspection accessors may be used to examine those results, but internal canonicalization, exceptional-frequency algorithms, and loop-coordinate machinery are not themselves part of the stable user API.

Coordinate-space Wick results are part of the qualified inspection surface through `Diagram`, `Diagrams`, `Edge`, `contractions`, `topology`, and `propagator_type`. Dressed propagators and self-energies expose their semantic R/A/K pieces through `retarded_component`, `advanced_component`, and `keldysh_component`; callers should not depend on representation storage fields.

Equal-time/Trotter provenance is controlled at the high-level propagator boundary with `DressedPropagator(...; preserve_regularisation=true)`. The underscored Wick-contraction control remains an implementation detail and is not part of the supported API.

The pre-rewrite `CollisionIntegral` path is legacy implementation code and is not part of the supported public API.
