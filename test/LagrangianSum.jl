using KeldyshContraction, Test
using KeldyshContraction: In, Out, Diagram, Diagrams, Edge
using KeldyshContraction: is_physical, is_conserved, _wick_contraction
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus
import KeldyshContraction as KC
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

@testset "Conversion" begin
    @test typeof(L) == KeldyshContraction.LagrangianSum{KeldyshContraction.QAdd{ComplexF64}}
end

@testset "Accessing" begin
    @test isequal(parameters(L), [Γ, g])
    # @test isequal(first(arguments(L)), L_inelastic)
end
