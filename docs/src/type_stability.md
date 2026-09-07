# Type-stability guarantee

KeldyshContraction.jl continuously enforces type stability for its explicitly supported
public API. The supported universe is deliberately bounded: a statistics/coefficient/API
combination is part of the guarantee only when it is represented in the public contract
workload and all compiler gates pass.

## Supported coefficient domains

For both supported statistics families, the guaranteed coefficient representations are:

- exact coefficients represented by `ComplexRationals`;
- floating-point coefficients represented by `ComplexF64`.

Other `Number` subtypes may work through generic methods, but they are outside the guarantee
until they are added to the supported-type contract matrix and CI workload.

## Bosonic public contract

For `Boson`, the guarantee covers package-owned public entry points for:

- field construction and access;
- coefficient conversion and rationalization;
- `InteractionLagrangian`;
- `wick_contraction`;
- `DressedPropagator` and `SelfEnergy`;
- matrix and topology access;
- Wigner transformation;
- package-owned plain-text and LaTeX rendering.

The resulting symbolic, diagram, propagator, self-energy, and Wigner representations must
remain recursively concrete and inference-visible.

## Fermionic public contract

For `Fermion`, the guarantee covers the v1 fermionic physics implemented by the LO layer:

- `FieldFamily{Fermion}` / `Field{Fermion}` construction through `@qfields`;
- `One` / `Two` component selection and `bar`;
- Grassmann field algebra;
- explicit coefficient conversion and rationalization;
- `InteractionLagrangian`, including multi-family target selection;
- `wick_contraction` with fermionic permutation signs;
- `DressedPropagator` and `SelfEnergy`;
- statistics-dispatched `matrix` access with the fermionic upper-triangular LO layout.

Representative exact and floating-point fermionic pipelines are exercised with `@inferred`
and recursively concrete result checks from field algebra through self-energy construction.

Fermionic Wigner transformation and collision-integral physics are not part of this guarantee.
They are explicitly outside the fermionic v1 scope and must not be inferred from the bosonic
coverage above.

## Enforcement

The guarantee is enforced by several independent gates:

1. the complete test suite runs with DispatchDoctor in `error` mode, so package-owned
   instability reached from tested supported workloads is an error;
2. every guaranteed coefficient/statistics domain is exercised through package-owned public
   computational APIs with `@inferred`;
3. resulting public representations are recursively checked for concrete field types;
4. JET package analysis must report no package-owned inference/runtime-dispatch errors;
5. `@unstable` exemptions are forbidden under `src/`;
6. Codecov is retained as an informational coverage trend and must upload successfully, but
   private implementation lines are not a type-stability acceptance target.

The coverage policy is deliberately API-oriented: tests are added to validate public behavior
and supported public type combinations, not merely to execute private branches for a
percentage. Internal tests may still validate implementation invariants when those invariants
are part of a physics or algebra correctness contract, but they do not define the public
coverage guarantee.

This is an engineering guarantee over the explicitly supported public input-type universe.
Extending it to another coefficient type, statistics family, physical transformation, or
public entry point requires adding that combination to the contract workload and making every
gate above pass.
