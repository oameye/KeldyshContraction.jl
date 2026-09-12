using KeldyshContraction
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus

#

@qfields ϕ::Boson
c, q = ϕ[Classical], ϕ[Quantum]
@syms Γ g

inelastic_terms =
    im * (
        0.5 * bar(c) * bar(q) * (c(Minus) * c(Minus) + q(Minus) * q(Minus)) -
        0.5 * c(Plus) * q(Plus) * (bar(c) * bar(c) + bar(q) * bar(q)) +
        bar(c) * bar(q) * (c(Plus) * q(Plus) + c(Minus) * q(Minus))
    )
elastic_terms = -(0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2))

L_inelastic = InteractionLagrangian(inelastic_terms, Γ)
L_elastic = InteractionLagrangian(elastic_terms, g)

L = L_inelastic + L_elastic

#

GF1 = DressedPropagator(L, Val(1), Val(3))

#

GF1_elastic = GF1[g]

#

GF1_inelastic = GF1[Γ]

#

GF2 = DressedPropagator(L, Val(2), Val(5); simplify=true)
topo = topologies(arguments(GF2)[g * Γ].keldysh)

#

[key => GF2[g * Γ].keldysh.diagrams[key] for key in topo[[2]]]

#

[key => GF2[g * Γ].keldysh.diagrams[key] for key in topo[[3]]]

#

Σ2 = SelfEnergy(GF2)
Σ2[g * Γ].keldysh

#

topo = topologies(Σ2[g * Γ].retarded)

[key => Σ2[g * Γ].retarded.diagrams[key] for key in topo[[2]]]

#

topo = topologies(Σ2[g * Γ].keldysh)
[key => Σ2[g * Γ].keldysh.diagrams[key] for key in topo[[3]]]
