# KeldyshContraction.jl

[![docs](https://img.shields.io/badge/docs-online-blue.svg)](https://oameye.github.io/KeldyshContraction.jl/)
[![codecov](https://codecov.io/gh/oameye/KeldyshContraction.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/oameye/KeldyshContraction.jl)
[![Benchmarks](https://github.com/oameye/KeldyshContraction.jl/actions/workflows/Benchmarks.yaml/badge.svg?branch=main)](https://oameye.github.io/KeldyshContraction.jl/benchmarks/)

[![Code Style: Blue](https://img.shields.io/badge/blue%20style%20-%20blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![jet](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)

KeldyshContraction.jl symbolically derives perturbative two-point functions, 1PI self-energies, and quantum-kinetic collision kernels in Schwinger--Keldysh field theory. The same pipeline supports bosonic and fermionic fields, coherent interactions, Lindblad loss vertices, coordinate derivatives, exact momentum routing, Wigner transformation, and statistics-dependent occupation reduction.

The main workflow is

```text
interaction action
    → DressedPropagator
    → SelfEnergy
    → Fourier transform
    → Wigner transform
    → Kadanoff–Baym collision expression
    → frequency/occupation reduction
    → CollisionKernel
```

The documentation is organised as a compact [theory introduction](https://oameye.github.io/KeldyshContraction.jl/dev/theory/keldysh/), a task-oriented [manual](https://oameye.github.io/KeldyshContraction.jl/dev/manual/fields/), and complete bosonic and fermionic transport examples.

The field-theory conventions follow the symbolic Keldysh framework developed in [`SciPost Phys. Core 8, 014 (2025)`](https://doi.org/10.21468/SciPostPhysCore.8.1.014).
