using KeldyshContraction
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus

#

@qfields c::Destroy(Classical) q::Destroy(Quantum)
@syms Γ g

inelastic_terms =
    im * (
        0.5 * c' * q' * (c(Minus) * c(Minus) + q(Minus) * q(Minus)) -
        0.5 * c(Plus) * q(Plus) * (c' * c' + q' * q') +
        c' * q' * (c(Plus) * q(Plus) + c(Minus) * q(Minus))
    )
elastic_terms = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))

L_inelastic = InteractionLagrangian(inelastic_terms, Γ)
L_elastic = InteractionLagrangian(elastic_terms, g)

L = L_inelastic + L_elastic

#

GF1 = DressedPropagator(L, 1; _set_reg_to_zero=true)

#

GF1_elastic = arguments(GF1)[g]

#

GF1_inelastic = arguments(GF1)[Γ]

#

GF2 = DressedPropagator(L, 2; _set_reg_to_zero=true)
topo = topologies(arguments(GF2)[g * Γ].keldysh)

#

topo[[2]]

#

[key => arguments(GF2)[g * Γ].keldysh.diagrams[key] for key in topo[[3]]]
