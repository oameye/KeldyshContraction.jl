# Keldysh-Schwinger Field Theory Conventions

This document fixes the Keldysh conventions used by KeldyshContraction.jl. Bosons use the
classical/quantum Retarded-Advanced-Keldysh basis. Fermions use the asymmetric
Larkin-Ovchinnikov rotation and the `One`/`Two` labels described below.

## Bosonic convention

The system is described by an action on the forward $(+)$ and backward $(-)$ contours. For
a bosonic field $\psi$ the contour fields are rotated to classical and quantum components,

```math
\psi_c = \frac{\psi_+ + \psi_-}{\sqrt{2}},\qquad
\psi_q = \frac{\psi_+ - \psi_-}{\sqrt{2}}.
```

The corresponding barred fields use the same rotation. For the quadratic diffusion sector,

```math
S_{\mathrm{diff}}^{RAK}
=\int dt\,d^dx\left\{
\bar\psi_q\left[i\partial_t+D\nabla^2-V(\mathbf{x})\right]\psi_c
+\bar\psi_c\left[i\partial_t+D\nabla^2-V(\mathbf{x})\right]\psi_q
\right\}.
```

The bosonic Green-function matrix is

```math
\hat G_B(x_1,x_2)
=
\begin{pmatrix}
G^K(x_1,x_2) & G^R(x_1,x_2)\\
G^A(x_1,x_2) & 0
\end{pmatrix}
=-i
\begin{pmatrix}
\langle\phi_c(x_1)\bar\phi_c(x_2)\rangle &
\langle\phi_c(x_1)\bar\phi_q(x_2)\rangle\\
\langle\phi_q(x_1)\bar\phi_c(x_2)\rangle &
\langle\phi_q(x_1)\bar\phi_q(x_2)\rangle
\end{pmatrix}.
```

The package uses the corresponding bosonic self-energy placement

```math
\hat\Sigma_B=
\begin{pmatrix}
0 & \Sigma^A\\
\Sigma^R & \Sigma^K
\end{pmatrix}.
```

`Classical` and `Quantum` are semantic aliases over the same neutral two-valued index stored
inside `Field{Boson}`; they do not create different Julia field types.

## Fermionic Larkin-Ovchinnikov convention

Fermions use the asymmetric Larkin-Ovchinnikov rotation. The unbarred and barred Grassmann
variables are rotated differently:

```math
\psi_1=\frac{\psi_+ + \psi_-}{\sqrt{2}},\qquad
\psi_2=\frac{\psi_+ - \psi_-}{\sqrt{2}},
```

```math
\bar\psi_1=\frac{\bar\psi_+ - \bar\psi_-}{\sqrt{2}},\qquad
\bar\psi_2=\frac{\bar\psi_+ + \bar\psi_-}{\sqrt{2}}.
```

The public labels `One` and `Two` denote these two components while reusing the same neutral
stored Keldysh index as the bosonic representation. Therefore a fermionic field remains a
`Field{Fermion}` regardless of its LO component.

With row index on the unbarred field and column index on the barred field, the package uses

```math
\hat G_F=
\begin{pmatrix}
G^R & G^K\\
0 & G^A
\end{pmatrix}_{1,2}.
```

Equivalently,

```text
(One, One) -> Retarded
(One, Two) -> Keldysh
(Two, One) -> structural zero
(Two, Two) -> Advanced
```

The fermionic self-energy uses the same upper-triangular LO placement,

```math
\hat\Sigma_F=
\begin{pmatrix}
\Sigma^R & \Sigma^K\\
0 & \Sigma^A
\end{pmatrix}.
```

This placement is part of the convention, not a relabeling of the bosonic classical/quantum
matrix. In particular, triangular matrices are closed under the Dyson product
$G_0\Sigma G_0$, so the structural-zero lower-left entry remains zero.

### `bar(psi)` and Grassmann variables

For fermions, `bar(psi)` denotes the independent barred Grassmann path-integral variable. It
toggles field orientation while preserving the physical field family, internal indices, LO
component, position, and regularisation.

It is deliberately not identified with the operator adjoint $\psi^\dagger$. The package does
not define `adjoint(::Field{Fermion}) = bar(field)`. Contraction- and edge-level adjoints are
separate operations used to relate retarded and advanced propagators.

### Fermionic Wick signs

Canonical fermionic products obey Grassmann algebra,

```math
A B=-B A,\qquad A^2=0
```

for distinct and identical fermionic generators, respectively. Generator identity includes
the field family and all identity-defining metadata, so different species or internal indices
do not vanish spuriously.

For Wick contraction, the package pairs the canonical unbarred sequence against the reversed
barred sequence. The pairing weight carries the parity of that permutation. This permutation
parity is the fermionic sign convention used by the implementation. No additional
$(-1)$ factor is applied per closed fermion loop; doing so would double-count a sign already
contained in the explicit Wick permutation parity.

## Perturbative Green functions and Dyson equation

For either statistics, the two-point dressed Green function is generated perturbatively from
an interaction action $S_{\mathrm{int}}=\int d^dx\,dt\,\mathcal L_{\mathrm{int}}$. Schematically,

```math
\hat G
=\hat G_0
+\hat G_0\circ\hat\Sigma\circ\hat G_0
+\hat G_0\circ\hat\Sigma\circ\hat G_0\circ\hat\Sigma\circ\hat G_0
+\cdots
=\hat G_0+\hat G_0\circ\hat\Sigma\circ\hat G.
```

Here $\circ$ denotes space-time convolution and Keldysh-matrix multiplication. The equivalent
Dyson equation is

```math
[\hat G_0^{-1}-\hat\Sigma]\hat G=\mathbb{1}.
```

Statistics dispatch determines the matrix placement of the stored retarded, advanced, and
Keldysh components; the shared `DressedPropagator` and `SelfEnergy` result types themselves
store those semantic components without duplicating the statistics-specific matrix layout.
