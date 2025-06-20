```@meta
EditURL = "../../../examples/collision_integral.jl"
```

````@example collision_integral
using KeldyshContraction
````

````@example collision_integral
@qfields c::Destroy(Classical) q::Destroy(Quantum)
elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
L_int = InteractionLagrangian(elasctic2boson)
````

````@example collision_integral
GF = DressedPropagator(L_int, 2)
````

````@example collision_integral
Σ = SelfEnergy(GF, 2)
````

````@example collision_integral
Σk = wigner_transform(Σ)
````

````@example collision_integral
ci = KeldyshContraction.CollisionIntegral(Σk)
ci.terms[[3]]
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

