# Type-stability guarantee

KeldyshContraction.jl continuously enforces type stability for its explicitly supported bosonic computational universe.

The current guarantee covers:

- `Boson` field statistics;
- exact coefficient storage represented by `ComplexRationals`;
- floating-point coefficient storage represented by `ComplexF64`;
- the package-owned pipeline from `InteractionLagrangian` through `DressedPropagator`, `SelfEnergy`, Wigner transformation, and `CollisionIntegral`;
- package-owned text and LaTeX serialization for objects exercised by that pipeline.

Other `Number` subtypes may work through generic methods, but they are outside this guarantee until they are added to the supported-type contract matrix and CI workload.

The guarantee is enforced by several independent gates:

1. the complete test suite runs with DispatchDoctor in `error` mode;
2. every supported coefficient domain is exercised through the public computational pipeline with `@inferred`;
3. resulting stored representations are recursively checked for concrete field types;
4. JET package analysis must report no package-owned inference/runtime-dispatch errors;
5. `@unstable` exemptions are forbidden under `src/`;
6. Codecov project and patch coverage targets are both 100%.

This is an engineering guarantee over the explicitly supported input-type universe. Extending the guarantee to another coefficient/statistics family requires adding that family to the contract workload and making all of the gates above pass.
