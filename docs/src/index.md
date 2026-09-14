# KeldyshContraction.jl

KeldyshContraction.jl symbolically derives perturbative two-point functions, 1PI self-energies, and quantum-kinetic collision kernels in Schwinger--Keldysh field theory. The same compiler supports bosonic and fermionic fields, coherent interactions, Lindblad loss vertices, coordinate derivatives, exact momentum routing, Wigner transformation, and statistics-dependent occupation reduction.

```@docs
KeldyshContraction
```

## Installation

The package is currently under active development. Install the development version directly from GitHub:

```julia
import Pkg
Pkg.add(url="https://github.com/oameye/KeldyshContraction.jl")
```

For a selected perturbative dressed-propagator sector `G`, the ordinary kinetic workflow is

```julia
C = collision_kernel(G)
```

The call runs the certified Fourier, self-energy, Wigner, Kadanoff--Baym, frequency, occupation, and momentum-reduction stages and returns the final strict-quasiparticle `CollisionKernel`. If a sector contains a genuine causal-frequency pinch or unresolved finite-Trotter state, the high-level compiler reports it instead of silently discarding it.

The individual compiler stages remain public for derivations and diagnostics. The canonical bosonic and fermionic examples deliberately unpack `collision_kernel(G)` one transformation at a time; their intermediate values render as physics-facing equations rather than internal storage dumps.

Start with the [theory](theory/keldysh.md) for conventions and approximations, the [manual](manual/fields.md) for the public workflow, or the complete bosonic and fermionic transport examples.

The field-theory conventions follow the symbolic Keldysh framework developed in [SciPost Phys. Core **8**, 014 (2025)](https://doi.org/10.21468/SciPostPhysCore.8.1.014). See [References](literature.md) for the supporting nonequilibrium and open-system literature.
