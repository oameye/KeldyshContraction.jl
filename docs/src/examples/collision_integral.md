```@meta
EditURL = "../../../examples/collision_integral.jl"
```

````@example collision_integral
using KeldyshContraction
````

````@example collision_integral
@qfields c::Boson(Classical) q::Boson(Quantum)
elasctic2boson = -(
    0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2)
)
L_int = InteractionLagrangian(elasctic2boson)
````

````@example collision_integral
GF = DressedPropagator(L_int, Val(2), Val(5))
````

````@example collision_integral
Σ = SelfEnergy(GF, Val(2))
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

