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

@testset "Correctness first order" begin
    GF1 = DressedPropagator(L, 1; simplify=false)
    GF1_elastic = arguments(GF1)[g]
    GF1_inelastic = arguments(GF1)[Γ]

    trued_elastic = DressedPropagator(L_elastic, 1; simplify=false)
    @test isequal(trued_elastic.keldysh, GF1_elastic.keldysh)
    @test isequal(trued_elastic.retarded, GF1_elastic.retarded)
    @test isequal(trued_elastic.advanced, GF1_elastic.advanced)

    trued_inelastic = DressedPropagator(L_inelastic, 1; simplify=false)
    @test isequal(trued_inelastic.keldysh, GF1_inelastic.keldysh)
    @test isequal(trued_inelastic.retarded, GF1_inelastic.retarded)
    @test isequal(trued_inelastic.advanced, GF1_inelastic.advanced)
end

@testset "Correctness second order" begin
    GF2 = DressedPropagator(L, 2; simplify=false)
    GF2_elastic = arguments(GF2)[g ^ 2]
    GF2_inelastic = arguments(GF2)[Γ ^ 2]

    trued_elastic = DressedPropagator(L_elastic, 2; simplify=false)
    @test isequal(trued_elastic.keldysh, GF2_elastic.keldysh)
    @test isequal(trued_elastic.retarded, GF2_elastic.retarded)
    @test isequal(trued_elastic.advanced, GF2_elastic.advanced)

    trued_inelastic = DressedPropagator(L_inelastic, 2; simplify=false)
    @test isequal(trued_inelastic.keldysh, GF2_inelastic.keldysh)
    @test isequal(trued_inelastic.retarded, GF2_inelastic.retarded)
    @test isequal(trued_inelastic.advanced, GF2_inelastic.advanced)

    @testset "cross-terms are swappable" begin
        term12 = L_elastic(1).lagrangian * L_inelastic(2).lagrangian
        term21 = L_elastic(2).lagrangian * L_inelastic(1).lagrangian
        regularise = KeldyshContraction.should_regularise(term12)
        diagrams12 = Diagrams{5,1}()
        diagrams21 = Diagrams{5,1}()
        for arg in arguments(term12)
            KeldyshContraction.wick_contraction!(
                diagrams12, c(Out()) * c'(In()) * arg; simplify=false, regularise
            )
        end
        for arg in arguments(term21)
            KeldyshContraction.wick_contraction!(
                diagrams21, c(Out()) * c'(In()) * arg; simplify=false, regularise
            )
        end
        @test isequal(diagrams21, diagrams12)
    end
end

@testset "SelfEnergy" begin
    GF2 = DressedPropagator(L, 2)
    Σ2 = SelfEnergy(GF2)
    # arguments(GF2)[g*Γ] |> typeof
end
