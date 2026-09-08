using KeldyshContraction
import KeldyshContraction as KC

@qfields probe_ψ::Fermion

ψ₁ = probe_ψ[One]
ψ₂ = probe_ψ[Two]
∂xψ₂ = partial(ψ₂, :x)
vertex = ψ₁ * ∂xψ₂ * bar(ψ₁) * bar(∂xψ₂)

L = InteractionLagrangian(vertex, :γ)
G = DressedPropagator(L, Val(1), Val(3); simplify=false)
Σ = SelfEnergy(G)

println("PWAVE_PROBE_BEGIN")
for (component_name, diagrams) in ((:retarded, Σ.retarded), (:keldysh, Σ.keldysh), (:advanced, Σ.advanced))
    println("COMPONENT=", component_name, " COUNT=", length(diagrams))
    for (diagram_index, pair) in enumerate(diagrams)
        diagram, coefficient = pair
        println("DIAGRAM=", diagram_index, " COEFF=", repr(coefficient))
        for (edge_index, edge) in enumerate(KC.contractions(diagram))
            out_field, in_field = KC.fields(edge)
            println(
                "EDGE=", edge_index,
                " POS=", repr(KC.positions(edge)),
                " TYPE=", KC.propagator_type(edge),
                " OUTDER=", repr(derivatives(out_field)),
                " INDER=", repr(derivatives(in_field)),
            )
        end

        routed = KC.FourierDiagram(diagram)
        println("ROUTING=", repr([m.coefficients for m in KC.edge_momenta(routed)]))
        lowered = KC.lower_fourier_derivatives(routed)
        polynomial = KC.kinematic_factor(lowered)
        println("POLY_TERMS=", length(polynomial.terms))
        for (monomial, value) in polynomial.terms
            factors = [(factor.axis, factor.momentum.coefficients) for factor in monomial.factors]
            println("TERM COEFF=", repr(value), " FACTORS=", repr(factors))
        end
    end
end
println("PWAVE_PROBE_END")
