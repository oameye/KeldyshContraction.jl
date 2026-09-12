# Supported workflow surface

The package-level public surface is intentionally smaller than the implementation IR.

The ordinary user workflow is exported and available directly after `using KeldyshContraction`. Representation and inspection utilities that are stable but too specialized for the caller namespace are marked public with SciMLPublic and are accessed as `KeldyshContraction.name` (or imported explicitly).

Internal frequency-reduction algorithms, canonicalization helpers, loop-coordinate transforms, and the legacy `CollisionIntegral` implementation are not API commitments.
