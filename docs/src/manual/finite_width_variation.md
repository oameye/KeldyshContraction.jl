# Finite-width spectral response

A finite-width spectral model may itself depend on the kinetic state. The first variation is therefore represented line by line as

```math
(E,\Gamma,Z)\longrightarrow(\delta E,\delta\Gamma,\delta Z).
```

For one Lorentzian line of multiplicity `m`,

```math
B_m(Z,\Gamma)=\binom{2m-2}{m-1}\frac{Z^m}{\Gamma^{m-1}},
```

and the package differentiates this exact expression directly. In particular,

```math
\delta\!\left(\frac{2Z^2}{\Gamma}\right)
=\frac{4Z}{\Gamma}\,\delta Z-
\frac{2Z^2}{\Gamma^2}\,\delta\Gamma.
```

Residual Cauchy convolutions reuse the already certified frequency geometry. Variations of line energies, linewidths, and residues are propagated through the exact mismatch, effective linewidth, and residue product. The microscopic `ω_external` is held fixed; it is not a collective center-time response frequency.

A supplied model variation is explicit. Missing required line variations are errors rather than implicit zeros.

```@docs
KeldyshContraction.LorentzianSpectralDataVariation
KeldyshContraction.LorentzianSpectralModelVariation
KeldyshContraction.lorentzian_integrated_power_variation
KeldyshContraction.evaluate_spectral_weight_variation
KeldyshContraction.evaluate_spectral_convolution_variation
```
