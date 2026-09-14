# Equal-time regularisation

Coincident-time propagators are not defined by the continuum contour notation alone. Their value depends on the operator ordering inherited from the discrete-time construction of the path integral.

For an equal-time contraction one must therefore retain whether a field approaches the time slice from above or below,

```math
t^+=t+0^+,
\qquad
t^-=t-0^+.
```

This information is physically relevant before the equal-time limit is taken.

## Trotter origin

In a coherent-state path integral, normally ordered operators act on neighbouring time slices after Trotter decomposition. Passing directly to a continuum action and setting all times equal erases that ordering information.

For closed-system vertices this often disappears after the appropriate cancellations. For Lindblad vertices it can change tadpoles, relative factors, and even whether disconnected contributions cancel. A dissipative action therefore has to retain the finite contour displacement implied by its microscopic discretisation.

KeldyshContraction.jl represents this provenance explicitly by `Plus`, `Zero`, and `Minus` regularisation labels on fields. These labels are carried through Wick contraction, Fourier transformation, and Wigner transformation.

## Relative equal-time shift

At the kinetic boundary the absolute endpoint labels are replaced by their physical relative shift,

```math
\delta r=r_{\mathrm{out}}-r_{\mathrm{in}}.
```

Thus propagators such as `G^R(y^+,y)` and `G^R(y,y^-)` represent the same regulated equal-time limit because both have `\delta r=+1`. This allows algebraically equivalent terms to combine without discarding their equal-time prescription.

The regulator is removed only when the corresponding structural equal-time reduction is performed. It is not absorbed into the generic definition of the spectral function or into an ad hoc numerical factor.

## Practical rule

Use

```julia
DressedPropagator(...; preserve_regularisation=true)
```

whenever the interaction contains finite-Trotter shifts, in particular for dissipative jump vertices. Ordinary interactions that do not depend on equal-time ordering can use the default.

The package does not invent the microscopic contour ordering of a jump operator. The interaction action must encode the correct finite-Trotter placement; the package then preserves and reduces that information consistently.

This distinction is essential in the two-body-loss examples, where the correct first-order loss kernel follows only after the regulated contour vertex is kept through the kinetic reduction.