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

## External spectral projection

The internal finite-width convolution is an off-shell collision weight at fixed microscopic external frequency. A momentum-occupation collision may require one further spectral projection over that external line. If

```math
W_{\rm int}(\omega)
=J Z_{\rm int}\frac{\Gamma_{\rm int}}
{(\beta\omega+\Delta_0)^2+\Gamma_{\rm int}^2/4},
```

then convolution with the normalized external Lorentzian is again exact:

```math
\int\frac{d\omega}{2\pi}A_{\rm ext}(\omega)W_{\rm int}(\omega)
=J Z_{\rm ext}Z_{\rm int}
\frac{\Gamma_{\rm full}}
{\Delta E^2+\Gamma_{\rm full}^2/4},
```

where

```math
\Gamma_{\rm full}=\Gamma_{\rm int}+|\beta|\Gamma_{\rm ext},
\qquad
\Delta E=\beta E_{\rm ext}+\Delta_0.
```

This is still a microscopic spectral integration. It is distinct from any collective center-time response frequency or moment-space pole condition. The exact directional derivative includes variations of the external line as well as every internal line.

```@docs
KeldyshContraction.external_spectral_projection_mismatch
KeldyshContraction.external_spectral_projection_linewidth
KeldyshContraction.evaluate_external_spectral_projection
KeldyshContraction.evaluate_external_spectral_projection_variation
```
