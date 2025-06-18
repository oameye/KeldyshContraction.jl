```@meta
EditURL = "../../../examples/two_particle_loss.jl"
```

## Two Body Loss

````@example two_particle_loss
using KeldyshContraction
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus
````

## System and regularization

The interaction action of elastic two body scattering, is defined as
```math
S_\mathrm{int} = i\Gamma \int d^d x \, [\frac{1}{2}(\bar{\phi}_+\phi_+)^2 +\frac{1}{2} (\bar{\phi}_-\phi_-)^2 -\bar{\phi}_-^2\phi_+^2]
```
Above interaction can typically represent s-wave scattering of bosons.

In the RAK basis, this gives
```math
S_\mathrm{int} = i\Gamma\int d^d x \, [\frac{1}{2}\bar{\phi}_c\bar{\phi}_q(\phi_c^2+\phi_q^2)-\frac{1}{2}\phi_c\phi_q(\bar{\phi}_c^2+\bar{\phi}_q^2)+2\bar{\phi}_c\bar{\phi}_q\phi_c\phi_q)]
```

A good check if the interaction Lagrangian is a valid physical process, is to check if the
normalization identity $Z=1$ holds. We can do this perturbatively in $g$ by expanding
$\exp(i S_\mathrm{int})$  and showing the average of the linear part of the system is zero
```math
\langle S_\mathrm{int}\rangle =  \langle S_\mathrm{int}^2\rangle  =\ldots = 0
```
In computing the average, one performs Wick contractions to describe the average in terms
of the two-point correlators of the linear part of the system. However, in this case we
don't find that the vacuum expectation value of the interaction Lagrangian is zero:

````@example two_particle_loss
@qfields c::Destroy(Classical) q::Destroy(Quantum)

loss2boson_unregular =
    im * (
        0.5 * c' * q' * (c^2 + q^2) - 0.5 * c * q * ((c')^2 + (q')^2) +
        c' * q' * (c * q + c * q)
    )

KeldyshContraction._wick_contraction(loss2boson_unregular; simplify=true)
````

To make the interaction Lagrangian physically meaningful, we must regularize it by properly
handling equal-space-time propagators. These equal-time arguments emerge from the continuum
limit of the field theory but are naturally absent in discrete-time formulations
[(Gerbino et al, 2024)](https://arxiv.org/abs/2406.20028).
In the discrete picture, operators act on coherent states at adjacent (but distinct) time
slices during path integral construction via Trotter decomposition. To ensure all disconnected
diagrams vanish identically, we introduce a finite time-shift regularization ε > 0 for the
quantum jump operators, motivated by the underlying discrete Trotter structure. One gets:
```math
\bar{L}_+(t)L_+(t) \to \bar{L}_+(t)L_+(t-\epsilon) \quad, \bar{L}_-(t)L_-(t) \to \bar{L}_-(t)L_-(t+\epsilon)
\quad \text{and} \quad
\bar{L}_-(t)L_+(t) \to \frac{1}{2}(\bar{L}_-(t)L_+(t+\epsilon) + \bar{L}_-(t)L_+(t-\epsilon))
```
Applying this regularization to the interaction Lagrangian, we get:

````@example two_particle_loss
loss2boson =
    im * (
        0.5 * c' * q' * (c(Minus) * c(Minus) + q(Minus) * q(Minus)) -
        0.5 * c(Plus) * q(Plus) * (c' * c' + q' * q') +
        c' * q' * (c(Plus) * q(Plus) + c(Minus) * q(Minus))
    )
L_int = InteractionLagrangian(loss2boson)
````

Indeed, the vacuum expectation value of the interaction Lagrangian is now zero:

## First order Green's function

````@example two_particle_loss
GF = DressedPropagator(L_int)
````

````@example two_particle_loss
Σ = SelfEnergy(GF)
````

The following indeed corresponds with what is reported in [(Gerbino et al, 2024)](https://arxiv.org/abs/2406.20028).

## Second order Green's function

````@example two_particle_loss
GF = DressedPropagator(L_int, 2)
````

````@example two_particle_loss
Σ = SelfEnergy(GF, 2)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

