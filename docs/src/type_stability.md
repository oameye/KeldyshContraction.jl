# Type-stability guarantee

KeldyshContraction.jl continuously enforces type stability for its explicitly supported bosonic public API.

The current guarantee covers:

- `Boson` field statistics;
- exact coefficient storage represented by `ComplexRationals`;
- floating-point coefficient storage represented by `ComplexF64`;
- package-owned public entry points for field construction and access, coefficient conversion, `InteractionLagrangian`, `wick_contraction`, `DressedPropagator`, `SelfEnergy`, matrix access, topology access, Wigner transformation, and package-owned rendering.

Other `Number` subtypes may work through generic methods, but they are outside this guarantee until they are added to the supported-type contract matrix and CI workload.

The guarantee is enforced by several independent gates:

1. the complete test suite runs with DispatchDoctor in `error` mode, so package-owned instability reached from tested public API workloads is an error;
2. every supported coefficient domain is exercised through the package-owned public computational API with `@inferred`;
3. resulting public representations are recursively checked for concrete field types;
4. JET package analysis must report no package-owned inference/runtime-dispatch errors;
5. `@unstable` exemptions are forbidden under `src/`;
6. Codecov is retained as an informational coverage trend and must upload successfully, but private implementation lines are not a type-stability acceptance target.

The coverage policy is deliberately API-oriented: tests are added to validate public behavior and supported public type combinations, not merely to execute private branches for a percentage. Codecov should remain high and should not regress without reason, but its repository-wide line percentage is not itself the type-stability guarantee.

This is an engineering guarantee over the explicitly supported public input-type universe. Extending the guarantee to another coefficient/statistics family or public entry point requires adding it to the contract workload and making all type-stability gates above pass.
