using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields canonical_collision_φ::Boson

function canonical_collision_test_term(
    kind::SpectralDispersiveKind, weight::StatisticalWeight; shift::Int8=Int8(0)
)
    basis = KC.MomentumBasis(2)
    q = KC.basis_momentum(basis, 2)
    line = KineticLine{Boson}(
        canonical_collision_φ, KC.Bulk(1), KC.Bulk(2), shift, q, KineticSpectral, weight
    )
    carrier = KineticTerm(
        KineticMonomial(KineticLine{Boson}[line], Val(1)),
        KC.FixedVector{0,Int}(Int[]),
        basis,
        basis[1],
        KC._kinematic_identity(),
    )
    return SpectralDispersiveTerm(carrier, KC.FixedVector{1,SpectralDispersiveKind}([kind]))
end

@testset "canonical collision frequency assembly" begin
    basis = KC.MomentumBasis(2)
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)

    Fk = KC.StatisticalAtom(canonical_collision_φ, k)
    Fq = KC.StatisticalAtom(canonical_collision_φ, q)
    @test KC.StatisticalMonomial([Fq, Fk]) == KC.StatisticalMonomial([Fk, Fq])

    spectral = canonical_collision_test_term(CollisionSpectral, DistributionWeight)
    spectral_expression = @inferred KC._causal_frequency_expression(
        spectral, 1 // 1, canonical_collision_φ
    )
    spectral_terms = KC.causal_frequency_terms(spectral_expression)
    @test length(spectral_terms) == 2
    @test Set(KC.causal_frequency_coefficient(term) for term in spectral_terms) ==
        Set((KC.ComplexRationals(im), -KC.ComplexRationals(im)))
    @test all(term -> length(KC.causal_frequency_denominators(term)) == 1, spectral_terms)

    dispersive = canonical_collision_test_term(CollisionDispersive, NoStatisticalWeight)
    dispersive_expression = @inferred KC._causal_frequency_expression(
        dispersive, 1 // 1, canonical_collision_φ
    )
    dispersive_terms = KC.causal_frequency_terms(dispersive_expression)
    @test length(dispersive_terms) == 2
    @test all(term -> KC.causal_frequency_coefficient(term) == 1 // 2, dispersive_terms)

    internal = KC._collision_statistical_monomial(spectral, canonical_collision_φ, false)
    complete = KC._collision_statistical_monomial(spectral, canonical_collision_φ, true)
    @test collect(internal) == [Fq]
    @test collect(complete) == sort([Fk, Fq])

    @testset "finite Trotter terms are preserved for the regulator layer" begin
        context = KC.HomogeneousWignerContext()
        shifted = canonical_collision_test_term(
            CollisionDispersive, NoStatisticalWeight; shift=Int8(1)
        )
        offset = SpectralDispersiveExpression{
            KC.ComplexRationals,Boson,1,0,0,KC.HomogeneousWignerContext
        }(
            context
        )
        distribution = SpectralDispersiveExpression{
            KC.ComplexRationals,Boson,1,0,0,KC.HomogeneousWignerContext
        }(
            context
        )
        push!(offset, shifted, 1 // 1)
        collision = SpectralDispersiveCollision{
            KC.ComplexRationals,Boson,1,1,0,0,KC.HomogeneousWignerContext
        }(
            offset, distribution, canonical_collision_φ, KC.ParameterMonomial(:γ), context
        )

        assembled = @inferred KC.canonical_frequency_collision(collision)
        @test isempty(KC.canonical_frequency_expressions(assembled))
        deferred = KC.shifted_frequency_terms(assembled)
        @test length(deferred) == 1
        @test KC.frequency_source(only(deferred)) == shifted
        @test KC.source_coefficient(only(deferred)) == 1
    end
end
