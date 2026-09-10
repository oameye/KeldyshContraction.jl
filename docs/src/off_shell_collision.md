# Off-shell collision identity

The spectral/statistical IR of `KineticSelfEnergy` is still an exact Wigner-space object. The next stage forms the zeroth-gradient Kadanoff--Baym collision combination, but does **not** yet make a quasiparticle or mass-shell approximation.

With the package convention

```math
A_\Sigma = i\left(\Sigma^R-\Sigma^A\right),
```

the homogeneous collision identity is

```math
I_{\mathrm{coll}}
= i\Sigma^K-iF_{\mathrm{target}}(k)\left(\Sigma^R-\Sigma^A\right)
= i\Sigma^K-F_{\mathrm{target}}(k)A_\Sigma.
```

The implementation stores this expression affinely in the external statistical distribution,

```math
I_{\mathrm{coll}}=I_0+F_{\mathrm{target}}(k)I_1,
\qquad
I_0=i\Sigma^K,
\qquad
I_1=-A_\Sigma.
```

This is the meaning of [`OffShellCollisionExpression`](@ref). `collision_offset` returns `I₀`; `collision_distribution_coefficient` returns `I₁`. The external distribution is not inserted as a `KineticLine`, because it belongs to the external leg rather than to the internal propagator product. Each retained `KineticTerm` therefore keeps its original topology, routed internal momenta, and derivative-generated momentum polynomial.

The spectral factor is formed from the **complete** retarded and advanced self-energies. In particular, the supported path is

```text
Σᴿ, Σᴬ
   ↓
Σᴿ - Σᴬ
   ↓
AΣ = i(Σᴿ - Σᴬ)
   ↓
iΣᴷ - Ftarget(k) AΣ
```

There is no per-edge replacement of the form `Im(product) -> product(Im(factor))`. For a genuine multi-line self-energy the complete `Σᴿ-Σᴬ` expression generally retains causal R/A products until the later spectral/dispersive decomposition.

## Greater/lesser cross-check

If greater and lesser self-energies are introduced with the standard algebraic relations

```math
\Sigma^K=\Sigma^>+\Sigma^<,
\qquad
\Sigma^R-\Sigma^A=\Sigma^>-\Sigma^<,
```

then the same collision identity becomes

```math
I_{\mathrm{coll}}
=i\left[(1-F)\Sigma^>+(1+F)\Sigma^<\right].
```

For bosons, `F = 1 + 2n`, giving

```math
I_{\mathrm{coll}}
=2i\left[(1+n)\Sigma^<-n\Sigma^>\right].
```

For fermions, `F = 1 - 2n`, giving

```math
I_{\mathrm{coll}}
=2i\left[n\Sigma^>+(1-n)\Sigma^<\right].
```

The statistics dependence therefore enters through the statistical convention rather than through a second collision-reduction algorithm.

## Target provenance

The physical external family is selected at `DressedPropagator(...; target=...)` for a multi-family theory and is preserved through self-energy amputation, Fourier transformation, Wigner transformation, and kinetic lowering. [`off_shell_collision_expression`](@ref) therefore takes only the `KineticSelfEnergy`; `target_family(KΣ)` is the authoritative external family used for `F_target(k)`.

This is important for same-statistics multi-family theories and for local self-energy sectors where the external family cannot be inferred reliably from the remaining internal lines.

## Approximation boundary

No operation on this page performs any of the following:

- replacing a spectral function by a shell delta function;
- imposing an energy-conservation constraint;
- replacing `F` by an occupation polynomial;
- resolving principal-value parts;
- integrating over phase space.

Those operations belong to the subsequent quasiparticle/on-shell projection. Keeping them separate is required so that the mixed `gγ` sector can retain its principal-value contribution instead of being forced into a delta-shell representation.

```@docs
OffShellCollisionExpression
off_shell_collision_expression
collision_offset
collision_distribution_coefficient
```
