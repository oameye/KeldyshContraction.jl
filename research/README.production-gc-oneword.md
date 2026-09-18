# Production GC one-word acceptance checkpoint

This branch promotes the certified `n <= 64` simple colored-directed GraphCombinations kernel into the physical dummy-loop quotient without runtime research patching.

Acceptance is defined by `.github/workflows/ProductionGCOneWordAcceptance.yml`:

- pin the exact GraphCombinations #204 production candidate;
- exercise `KC.quotient_loop_momenta` directly;
- compare against the independent Nauty projective-loop oracle;
- certify all 160 two-loop and all 160 three-loop production-shaped cases;
- require zero semantic failures and zero memory regressions;
- confirm any timing regression with the existing repeated benchmark protocol.

The production path uses the public GC `DirectedSimpleCanonicalizationWorkspace` / `canonicalize_directed_simple!` boundary and applies the resulting signed loop permutation directly to occupations, kinematic momenta, and shell/PV support. Nauty remains an independent oracle, not the production backend.
