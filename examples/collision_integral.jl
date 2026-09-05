using KeldyshContraction

#

@qfields ϕ::Boson
c, q = ϕ[Classical], ϕ[Quantum]
elasctic2boson = -(
    0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2)
)
L_int = InteractionLagrangian(elasctic2boson)

#

GF = DressedPropagator(L_int, Val(2), Val(5))

#

Σ = SelfEnergy(GF, Val(2))

#

Σk = wigner_transform(Σ)

#

ci = KeldyshContraction.CollisionIntegral(Σk)
ci.terms[[3]]
