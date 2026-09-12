using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields derivative_ϕ::Boson derivative_ψ::Fermion

@testset "canonical derivative metadata" begin
    c = derivative_ϕ[Classical]
    ψ₁ = derivative_ψ[One]

    ∂xψ = @inferred partial(ψ₁, :x)
    ∂yψ = @inferred partial(ψ₁, :y)
    ∂xyψ = @inferred partial(partial(ψ₁, :x), :y)
    ∂yxψ = @inferred partial(partial(ψ₁, :y), :x)
    ∂xxψ = @inferred partial(partial(ψ₁, :x), :x)
    ∂tzψ = @inferred partial(partial(ψ₁, :z), :t)

    @test isbitstype(KC.DerivativeMultiIndex)
    @test @inferred(derivatives(ψ₁)) == Symbol[]
    @test @inferred(derivatives(∂xψ)) == [:x]
    @test @inferred(derivatives(∂yψ)) == [:y]
    @test @inferred(derivatives(∂xyψ)) == [:x, :y]
    @test @inferred(derivatives(∂xxψ)) == [:x, :x]
    @test @inferred(derivatives(∂tzψ)) == [:t, :z]
    @test isequal(∂xyψ, ∂yxψ)
    @test hash(∂xyψ) == hash(∂yxψ)
    @test_throws ArgumentError partial(ψ₁, :r)

    # Public derivative access is owned: mutating the returned vector cannot mutate the field.
    axes = derivatives(∂xyψ)
    push!(axes, :z)
    @test derivatives(∂xyψ) == [:x, :y]

    @test field_family(∂xψ) === derivative_ψ
    @test field_family(partial(c, :x)) === derivative_ϕ
end

@testset "statistics-neutral derivative algebra" begin
    c = derivative_ϕ[Classical]
    ψ₁ = derivative_ψ[One]
    ∂xc = @inferred partial(c, :x)
    ∂xψ = @inferred partial(ψ₁, :x)

    @test !isequal(c, ∂xc)
    @test isequal(@inferred(c * ∂xc), @inferred(∂xc * c))

    @test !iszero(@inferred(ψ₁ * ∂xψ))
    @test isequal(@inferred(ψ₁ * ∂xψ), @inferred(-1 * (∂xψ * ψ₁)))
    @test iszero(@inferred(∂xψ * ∂xψ))

    ∂xxψ = @inferred partial(∂xψ, :x)
    @test iszero(@inferred(∂xxψ * ∂xxψ))
    @test !iszero(@inferred(∂xψ * ∂xxψ))
end

@testset "field reconstruction preserves derivatives" begin
    using KeldyshContraction: Bulk, Regularisation

    ψ₁ = derivative_ψ[One]
    ∂xψ = partial(ψ₁, :x)

    @test derivatives(@inferred(bar(∂xψ))) == [:x]
    @test derivatives(@inferred(∂xψ(Bulk(3)))) == [:x]
    @test derivatives(@inferred(∂xψ(Regularisation.Plus))) == [:x]
    @test derivatives(@inferred(KC.set_reg_to_zero(∂xψ(Regularisation.Plus)))) == [:x]
end

@testset "derivatives are passive for contraction physics" begin
    using KeldyshContraction: Contraction, PropagatorType

    ψ₁ = derivative_ψ[One]
    ψ₂ = derivative_ψ[Two]
    ∂xψ₁ = partial(ψ₁, :x)
    ∂xψ₂ = partial(ψ₂, :x)

    @test KC.contraction_compatible(Fermion, ∂xψ₁, bar(ψ₁))
    @test KC.contraction_compatible(Fermion, ψ₁, bar(∂xψ₁))
    @test KC.propagator_type(Fermion, ∂xψ₁, bar(ψ₁)) === PropagatorType.Retarded
    @test KC.propagator_type(Fermion, ∂xψ₁, bar(∂xψ₂)) === PropagatorType.Keldysh

    left = Contraction(∂xψ₁, bar(ψ₁))
    right = Contraction(ψ₁, bar(∂xψ₁))
    @test !isequal(KC.propagator_color(left), KC.propagator_color(right))

    canonical_left = only(@inferred KC.canonicalize([left]))
    canonical_right = only(@inferred KC.canonicalize([right]))
    @test !isequal(canonical_left, canonical_right)
end

@testset "spinless derivative interaction reaches self energy" begin
    ψ₁ = derivative_ψ[One]
    ψ₂ = derivative_ψ[Two]
    ∂xψ₂ = partial(ψ₂, :x)

    vertex = @inferred ψ₁ * ∂xψ₂ * bar(ψ₁) * bar(∂xψ₂)
    @test !iszero(vertex)

    L = @inferred InteractionLagrangian(vertex, :γ)
    G = @inferred DressedPropagator(L, Val(1), Val(3); simplify=false)
    Σ = @inferred SelfEnergy(G)

    @test G isa DressedPropagator{KC.ComplexRationals,Fermion,1,3,0}
    @test Σ isa SelfEnergy{KC.ComplexRationals,Fermion,1,1,0}

    function contains_derivative(ds)
        for diagram in keys(ds.diagrams),
            edge in KC.contractions(diagram),
            field in KC.fields(edge)

            isempty(derivatives(field)) || return true
        end
        return false
    end

    @test contains_derivative(G.retarded) ||
        contains_derivative(G.keldysh) ||
        contains_derivative(G.advanced)
    @test contains_derivative(Σ.retarded) ||
        contains_derivative(Σ.keldysh) ||
        contains_derivative(Σ.advanced)
end
