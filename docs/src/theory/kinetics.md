# Quantum kinetic theory

Quantum transport follows from the Keldysh Dyson equation after separating spectral information from the distribution function. KeldyshContraction.jl focuses on the collision side of this derivation and keeps each approximation explicit.

## From Dyson to the kinetic equation

Write the exact Keldysh propagator as

```math
G^K=G^R\circ F-F\circ G^A,
```

where `\circ` denotes space-time convolution. Combining the left and right Dyson equations gives

```math
[F,G_0^{-1}]_\circ
=\Sigma^K-\Sigma^R\circ F+F\circ\Sigma^A.
```

The left-hand side generates streaming and gradient corrections. The right-hand side is the collision term. This identity is off shell and does not assume a quasiparticle form.

## Wigner representation

For a two-point function `G(x_1,x_2)`, introduce centre and relative coordinates,

```math
X=\frac{x_1+x_2}{2},
\qquad
s=x_1-x_2,
```

and Fourier transform in `s`. Convolutions become Moyal products. Expanding the Moyal product converts the left-hand commutator into the familiar drift/Poisson-bracket structure. The current collision compiler uses the homogeneous zeroth-gradient limit for the right-hand side, where the products reduce to ordinary multiplication.

With

```math
A_\Sigma=i\left(\Sigma^R-\Sigma^A\right),
```

the collision combination becomes

```math
I_{\mathrm{coll}}
=i\Sigma^K-iF_k\left(\Sigma^R-\Sigma^A\right)
=i\Sigma^K-F_kA_\Sigma.
```

It is formed from the **complete** retarded and advanced self-energies. Products of causal propagators are therefore reduced only after the full Kadanoff--Baym combination has been assembled.

For internal lines at zeroth gradient order,

```math
G^K=-iFA,
```

with the statistics-dependent occupation maps

```math
F_B=1+2n_B,
\qquad
F_F=1-2n_F.
```

The collision-reduction algorithm itself is consequently shared by bosons and fermions; statistics enters through the field algebra, Wick signs, and the final map from `F` to `n`.

## Spectral and dispersive parts

A causal propagator carries both spectral and dispersive information. Frequency reduction therefore retains both

```math
\delta(\Delta E),
\qquad
\operatorname{PV}\!\left(\frac{1}{\Delta E}\right).
```

The shell terms generate on-shell gain/loss processes. Principal-value terms describe virtual dispersive contributions and are not discarded merely because they are off shell. Momentum conservation is exact; equivalent loop-coordinate choices are quotiented only after the occupation polynomial is formed.

## Strict quasiparticle reduction

In the strict quasiparticle limit the spectral support collapses onto single-particle shells. After the external frequency is projected onto its physical shell, the finite part takes the schematic form

```math
C[n](k)
=\int_{\boldsymbol q\ldots}
\mathcal K(k,\boldsymbol q,\ldots)
\,\mathcal N[n]
\,\delta(\Delta E),
```

with principal-value sectors retained when present. `\mathcal K` is the exact momentum polynomial produced by the vertices and `\mathcal N[n]` is the bosonic or fermionic occupation polynomial.

Not every causal structure has a finite strict-quasiparticle limit. Degenerate retarded/advanced poles can produce a genuine pinch singularity. Such terms remain explicit blockers; KeldyshContraction.jl does **not** assign them an arbitrary finite coefficient. Finite-width or resummed kinetics is a separate approximation.

### Bosonic two-body-loss second-order oracle

For the local bosonic jump operator `L = ψ²`, the analytical second-order benchmark used by the test suite is the two-part result labelled Eq. (55a,b) in the reference derivation:

```math
I_{\mathrm{coll}}^{(\gamma^2)}[n]
=
4\gamma^2\int_{q_1,q_2,q_3}
\mathcal N_{55a}[n]\,
A_{q_1}A_{q_2}A_{q_3}
+
16\gamma^2\int_{q_1,q_2}
n_k n_{q_1}n_{q_2}\,
A_{q_1}^2A_{q_2},
```

where

```math
\mathcal N_{55a}[n]
=
(1+n_k)(1+n_{q_3})n_{q_1}n_{q_2}
+(1+n_{q_1})(1+n_{q_2})n_{q_3}n_k
-2n_{q_1}n_{q_2}n_{q_3}n_k.
```

The package stores the occupation-number equation rather than the equation for `F_B = 1 + 2n_B`, so for this bosonic convention

```math
C_n=\frac{1}{2}I_{\mathrm{coll}}.
```

Consequently Eq. (55a) becomes the regular `2γ² 𝒩₅₅ₐ[n]` kernel. Expanding `𝒩₅₅ₐ` cancels the quartic occupation term and gives the six-term polynomial checked by `generated_collision_kernel.jl`.

Eq. (55b) is the second piece of the **same** `O(γ²)` result. In the package normalization it carries an `8γ²` prefactor, but its repeated spectral line `A_{q_1}²` is a degenerate retarded/advanced shell. The strict-quasiparticle compiler therefore certifies this contribution as repeated-shell pinch support in `generated_loss_canonical_support.jl` and deliberately keeps it upstream of `CollisionKernel` rather than assigning it a spurious finite coefficient.

The package derives `C[n]`; it does not currently construct or solve the streaming equation on the left-hand side.

## Compiler chain

The ordinary interface hides the mechanics without hiding the approximation:

```text
DressedPropagator sector G
  ↓ collision_kernel
C[n]
```

For derivations and diagnostics, the same call can be unpacked into the certified compiler:

```text
G
  ↓ fourier_transform
Fourier propagator
  ↓ SelfEnergy
Σᴿ,Σᴬ,Σᴷ
  ↓ wigner_transform
Wigner self-energy
  ↓ kinetic_expression
spectral/statistical self-energy
  ↓ off_shell_collision_expression
Kadanoff--Baym collision expression
  ↓ spectral_dispersive_collision
spectral/dispersive basis
  ↓ reduce_frequency_collision
shell / PV / blocked frequency support
  ↓ occupation_reduced_expression
occupation polynomial
  ↓ quotient_loop_momenta
canonical loop representation
  ↓ collision_kernel
C[n]
```

The [Kinetic reduction](../manual/kinetics.md) manual page documents the corresponding public functions, while the canonical examples execute each stage and render its physics directly.
