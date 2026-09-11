using KeldyshContraction
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus

#

@qfields c::Boson(Classical) q::Boson(Quantum)
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

GF1 = DressedPropagator(L, Val(1), Val(3); _set_reg_to_zero=true)

#

GF1_elastic = arguments(GF1)[g]

#

GF1_inelastic = arguments(GF1)[Γ]

#

GF2 = DressedPropagator(L, Val(2), Val(5); _set_reg_to_zero=true, simplify=true)
topo = topologies(arguments(GF2)[g * Γ].keldysh)

#

[key => arguments(GF2)[g * Γ].keldysh.diagrams[key] for key in topo[[2]]]

#

[key => arguments(GF2)[g * Γ].keldysh.diagrams[key] for key in topo[[3]]]

#

Σ2 = SelfEnergy(GF2, Val(2))
arguments(Σ2)[g * Γ].keldysh

#

topo = topologies(arguments(Σ2)[g * Γ].retarded)

[key => arguments(Σ2)[g * Γ].retarded.diagrams[key] for key in topo[[2]]]

#

topo = topologies(arguments(Σ2)[g * Γ].keldysh)
[key => arguments(Σ2)[g * Γ].keldysh.diagrams[key] for key in topo[[3]]]
