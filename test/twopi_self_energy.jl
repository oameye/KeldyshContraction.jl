using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields twopi_derivative_ϕ::Boson

function twopi_derivative_interaction()
    c = twopi_derivative_ϕ[Classical]
    q = twopi_derivative_ϕ[Quantum]
    return InteractionLagrangian(c * q * bar(c) * bar(q), :λ)
end

@testset "formal 2PI graph derivative precedes physical projection" begin
    L = twopi_derivative_interaction()
    Γ2 = TwoPIEffectiveAction(L, Val(1), Val(2))
    Σformal = @inferred KC.twopi_self_energy(Γ2, twopi_derivative_ϕ)

    @test KC.order(Σformal) == 1
    @test KC.statistics(Σformal) === Boson
    @test KC.parameters(Σformal) == KC.parameters(Γ2)
    @test KC.target_family(Σformal) == twopi_derivative_ϕ

    terms = KC.twopi_self_energy_terms(Σformal)
    @test length(terms) == 4
    @test all(value == -im for value in values(terms))

    keldysh_cuts = [
        diagram for diagram in keys(terms) if
        KC.twopi_self_energy_component(diagram) == (Quantum, Quantum)
    ]
    @test length(keldysh_cuts) == 1
    @test !any(KC.is_qq_contraction, KC.twopi_self_energy_contractions(only(keldysh_cuts)))

    structural_zero_cuts = [
        diagram for diagram in keys(terms) if
        KC.twopi_self_energy_component(diagram) == (Classical, Classical)
    ]
    @test length(structural_zero_cuts) == 1
    @test any(
        KC.is_qq_contraction, KC.twopi_self_energy_contractions(only(structural_zero_cuts))
    )

    # The physical Keldysh self-energy exists because the formal qq line was retained in Γ₂ and
    # cut before the physical Gqq = 0 restriction was applied.
    Σphysical = @inferred KC.physical_self_energy(Σformal)
    @test !iszero(KC.keldysh_component(Σphysical))
end

@testset "2PI derivative matches independent skeleton self-energy" begin
    L = twopi_derivative_interaction()
    Γ2 = TwoPIEffectiveAction(L, Val(1), Val(2))
    projected = KC.physical_self_energy(KC.twopi_self_energy(Γ2, twopi_derivative_ϕ))

    G = DressedPropagator(L, Val(1), Val(3); simplify=false)
    reference = KC.skeleton_self_energy(G)

    @test KC.parameters(projected) == KC.parameters(reference)
    @test KC.target_family(projected) == KC.target_family(reference)
    @test KC.keldysh_component(projected) == KC.keldysh_component(reference)
    @test KC.retarded_component(projected) == KC.retarded_component(reference)
    @test KC.advanced_component(projected) == KC.advanced_component(reference)
end

@qfields twopi_derivative_ψ::Fermion

@testset "fermionic derivative normalization is explicit" begin
    c = twopi_derivative_ψ[Classical]
    q = twopi_derivative_ψ[Quantum]
    L = InteractionLagrangian(c * q * bar(c) * bar(q), :u)
    Γ2 = TwoPIEffectiveAction(L, Val(1), Val(2))

    @test_throws ArgumentError KC.twopi_self_energy(Γ2, twopi_derivative_ψ)
end
