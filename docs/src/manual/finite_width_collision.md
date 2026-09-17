# Finite-width collision reduction

The finite-width branch acts on the complete spectral/dispersive collision before the sharp-quasiparticle frequency limit. It preserves the exact affine Kadanoff--Baym structure

```math
I_{\mathrm{coll}}=I_0+F_{\mathrm{target}}(k)I_1
```

while classifying each routed term with the analytic Lorentzian backends. Independent or repeated spectral lines use the factorized reduction, the regular `L+1`-line / `L`-loop sector uses the residual Cauchy convolution, and structures outside those analytic domains remain explicitly unsupported.

This stage does not supply a linewidth model, choose the microscopic external spectral frequency, lower `F` to occupations, or integrate momenta. In particular, the current statistical factors remain the quasiparticle/statistical factors carried by the kinetic IR; this is not a general solver for a frequency-dependent distribution `F(ω,k)`.

```@docs
KeldyshContraction.FiniteWidthFrequencyReductionKind
KeldyshContraction.FiniteWidthFrequencyTerm
KeldyshContraction.FiniteWidthFrequencyExpression
KeldyshContraction.FiniteWidthFrequencyCollision
KeldyshContraction.finite_width_frequency_collision
```
