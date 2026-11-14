using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields twopi_parity_ϕ::Boson

@testset "single-family 2PI route reaches ordinary SelfEnergy" begin
    c = twopi_parity_ϕ[Classical]
    q = twopi_parity_ϕ[Quantum]
    L = InteractionLagrangian(c * q * bar(c) * bar(q), :λ)

    Γ2 = TwoPIEffectiveAction(L, Val(1), Val(2))
    Σ2 = @inferred SelfEnergy(Γ2)

    G = DressedPropagator(L, Val(1), Val(3); simplify=false)
    reference = KC.skeleton_self_energy(G)

    @test Σ2 isa SelfEnergy
    @test KC.order(Σ2) == KC.order(reference)
    @test KC.statistics(Σ2) === KC.statistics(reference)
    @test KC.parameters(Σ2) == KC.parameters(reference)
    @test KC.target_family(Σ2) == KC.target_family(reference)
    @test KC.keldysh_component(Σ2) == KC.keldysh_component(reference)
    @test KC.retarded_component(Σ2) == KC.retarded_component(reference)
    @test KC.advanced_component(Σ2) == KC.advanced_component(reference)
    @test KC.matrix(Σ2) == KC.matrix(reference)
end

@qfields twopi_parity_ψ::Boson twopi_parity_χ::Boson

function twopi_parity_hs_interaction()
    ψc = twopi_parity_ψ[Classical]
    ψq = twopi_parity_ψ[Quantum]
    χq = twopi_parity_χ[Quantum]

    forward = ψc^2 * bar(χq) + 2 * ψc * ψq * bar(χq) + (1 // 4) * ψq^2 * bar(χq)
    interaction = -im * (forward + bar(forward))
    return ChargedInteractionLagrangian(
        interaction, twopi_parity_ψ => 1, twopi_parity_χ => 2; parameter=:h
    )
end

function test_same_self_energy_semantics(actual, expected)
    @test KC.order(actual) == KC.order(expected)
    @test KC.statistics(actual) === KC.statistics(expected)
    @test KC.parameters(actual) == KC.parameters(expected)
    @test KC.target_family(actual) == KC.target_family(expected)
    @test KC.keldysh_component(actual) == KC.keldysh_component(expected)
    @test KC.retarded_component(actual) == KC.retarded_component(expected)
    @test KC.advanced_component(actual) == KC.advanced_component(expected)
    @test KC.matrix(actual) == KC.matrix(expected)
    return nothing
end

@testset "multi-family HS 2PI route reaches ordinary SelfEnergy" begin
    L = twopi_parity_hs_interaction()
    Γ2 = TwoPIEffectiveAction(L, Val(2), Val(3))

    @test_throws ArgumentError SelfEnergy(Γ2)

    Σψ = @inferred SelfEnergy(Γ2, twopi_parity_ψ)
    Ωχ = @inferred SelfEnergy(Γ2, twopi_parity_χ)
    explicit_Σψ = @inferred KC.physical_self_energy(
        KC.twopi_self_energy(Γ2, twopi_parity_ψ)
    )
    explicit_Ωχ = @inferred KC.physical_self_energy(
        KC.twopi_self_energy(Γ2, twopi_parity_χ)
    )

    @test Σψ isa SelfEnergy
    @test Ωχ isa SelfEnergy
    test_same_self_energy_semantics(Σψ, explicit_Σψ)
    test_same_self_energy_semantics(Ωχ, explicit_Ωχ)

    @test KC.order(Σψ) == 2
    @test KC.order(Ωχ) == 2
    @test KC.statistics(Σψ) === Boson
    @test KC.statistics(Ωχ) === Boson
    @test KC.parameters(Σψ) == KC.parameters(Γ2)
    @test KC.parameters(Ωχ) == KC.parameters(Γ2)
    @test KC.target_family(Σψ) == twopi_parity_ψ
    @test KC.target_family(Ωχ) == twopi_parity_χ

    Mψ = @inferred KC.matrix(Σψ)
    Mχ = @inferred KC.matrix(Ωχ)
    @test size(Mψ) == (2, 2)
    @test size(Mχ) == (2, 2)
    @test Mψ[1, 2] == KC.advanced_component(Σψ)
    @test Mψ[2, 1] == KC.retarded_component(Σψ)
    @test Mψ[2, 2] == KC.keldysh_component(Σψ)
    @test Mχ[1, 2] == KC.advanced_component(Ωχ)
    @test Mχ[2, 1] == KC.retarded_component(Ωχ)
    @test Mχ[2, 2] == KC.keldysh_component(Ωχ)
end
