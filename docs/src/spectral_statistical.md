# Spectral and statistical kinetic representation

The kinetic layer starts from the semantic retarded, advanced, and Keldysh components of a `WignerSelfEnergy`. It does not infer physics from the matrix position of those components. This matters because the package uses different bosonic c/q and fermionic Larkin-Ovchinnikov matrix layouts while exposing the same semantic `retarded`, `advanced`, and `keldysh` fields.

## Spectral convention

The package kinetic convention defines the spectral function by

```math
A = i\left(G^R-G^A\right).
```

Equivalently,

```math
G^R-G^A=-iA.
```

The same normalization is used for a complete self-energy,

```math
A_\Sigma=i\left(\Sigma^R-\Sigma^A\right).
```

The subtraction is performed on the complete retarded and advanced self-energy expressions. It is not replaced by a product of imaginary parts of individual propagator lines.

## Statistical function

Before a gradient approximation, the exact kinetic parametrisation is

```math
G^K = G^R\circ F-F\circ G^A,
```

where `\circ` denotes the same convolution/Moyal product inherited from the Wigner representation. In the currently supported homogeneous gradient-order-zero regime, `F` commutes with the propagators at the algebraic level, so

```math
G^K=F\left(G^R-G^A\right)=-iFA.
```

This identity is the only per-line R/A/K basis conversion performed by `kinetic_expression`. Every Keldysh propagator becomes one spectral line carrying one statistical weight `F` and contributes one exact factor `-i` to the numerical coefficient. Retarded and advanced propagators remain retarded and advanced.

The previous experimental collision reducer used `G^K=(i/2)FA`. That normalization is not used by the new kinetic IR.

## Occupation convention

For bosons, with

```math
G^> = -i(1+n_B)A,
\qquad
G^< = -i\,n_B A,
```

one obtains

```math
G^K=G^>+G^<=-i(1+2n_B)A,
```

and therefore

```math
F_B=1+2n_B.
```

For fermions, with

```math
G^> = -i(1-n_F)A,
\qquad
G^< = +i\,n_F A,
```

one obtains

```math
G^K=G^>+G^<=-i(1-2n_F)A,
```

and therefore

```math
F_F=1-2n_F.
```

The shared kinetic algebra therefore needs statistics dispatch only for the affine occupation map. Spectral storage, exact momentum routing, kinematic polynomials, duplicate merging, and complete-component subtraction are statistics-neutral.

## Concrete IR

`KineticLine{S}` stores one physical propagator line after the Keldysh matrix index has been eliminated. It retains:

- the physical `FieldFamily{S}`;
- out/in positions;
- the canonical relative equal-time regularisation shift;
- the exact routed `LinearMomentum`;
- a causal/spectral `KineticLineKind`;
- a `StatisticalWeight` specifying whether the line carries `F`.

`KineticMonomial{S,E}` is the canonical commutative product of `E` kinetic lines. It has a statically visible line count and canonical line order, so products that become identical only after R/A/K-to-spectral/statistical conversion merge deterministically rather than inheriting an accidental source-diagram edge ordering.

`KineticTerm{S,E1,E2}` combines that monomial with the fixed topology signature, complete `MomentumBasis`, external Wigner momentum, and exact derivative-generated `MomentumPolynomial`.

`KineticExpression` is a concrete dictionary from canonical kinetic terms to numeric coefficients. Equal terms merge exactly and zero coefficients are removed. `KineticSelfEnergy` keeps the complete Keldysh, retarded, and advanced expressions separate.

The complete self-energy discontinuity is constructed as

```julia
KΣ = kinetic_expression(ΣW)
ΔΣ = retarded_minus_advanced(KΣ)
AΣ = spectral_self_energy(KΣ)
```

No quasiparticle projection, shell delta function, phase-space integration, or final collision expression is introduced at this stage.

## Equal-time regularisation is separate

The spectral definition above must not be confused with the special equal-time limit required by dissipative tadpole diagrams. Coordinate/Fourier/Wigner objects preserve the raw endpoint `Regularisation.Plus`, `Regularisation.Zero`, and `Regularisation.Minus` tags. At the kinetic boundary those raw spellings are replaced by the physical relative shift

```math
\delta r = r_{\rm out}-r_{\rm in}.
```

For example, `Gᴿ(y⁺,y)` and `Gᴿ(y,y⁻)` both have `δr=+1` and therefore become the same `KineticLine`. This is the correct semantic identity for the regulated equal-time propagator and allows duplicate tadpole contributions to merge. The sign and magnitude of the relative shift remain explicit for the later equal-time reduction.

The ordinary `DressedPropagator` path already has the required provenance boundary. For a regularised interaction, `_set_reg_to_zero=true` removes regulator decoration from different-position propagators and from Keldysh propagators, where it is not needed downstream, but retains the regulator on same-position retarded/advanced propagators. Fourier and Wigner preserve those tags; kinetic lowering canonicalizes them to `regularisation_shift`.

Any factor-of-two from the regulated equal-time limit belongs to the later reduction of those specific terms. It is not folded into the definition of `A` or into generic R/A/K conversion.

## Lesser/greater basis

The equations above fix lesser/greater conventions unambiguously, but the initial kinetic IR does not introduce a second lesser/greater diagram container. If a later gain/loss rewrite benefits from explicit lesser/greater components, it should be implemented as a semantic conversion over the same `KineticExpression`, with tests against the R/A/K identities above.
