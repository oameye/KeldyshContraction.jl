using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields skeleton_loss_ϕ::Boson
@qfields skeleton_pwave_ψ::Fermion

function skeleton_loss_lagrangian()
    c = skeleton_loss_ϕ[Classical]
    q = skeleton_loss_ϕ[Quantum]
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus
    loss =
        (1 // 2) *
        im *
        (
            bar(c) * bar(q) * (c(minus)^2 + q(minus)^2) -
            c(plus) * q(plus) * (bar(c)^2 + bar(q)^2) +
            2 * bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
        )
    return InteractionLagrangian(loss, :γ)
end

function skeleton_pwave_lagrangian()
    ψ₁ = skeleton_pwave_ψ[One]
    ψ₂ = skeleton_pwave_ψ[Two]
    bψ₁ = bar(ψ₁)
    bψ₂ = bar(ψ₂)
    ∂bψ₁ = partial(bψ₁, :x)
    ∂bψ₂ = partial(bψ₂, :x)
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus

    ψplus_minus = ψ₁(minus) + ψ₂(minus)
    ∂ψplus_minus = partial(ψ₁(minus), :x) + partial(ψ₂(minus), :x)
    ψplus_plus = ψ₁(plus) + ψ₂(plus)
    ∂ψplus_plus = partial(ψ₁(plus), :x) + partial(ψ₂(plus), :x)
    ψminus_plus = ψ₁(plus) - ψ₂(plus)
    ∂ψminus_plus = partial(ψ₁(plus), :x) - partial(ψ₂(plus), :x)

    bψplus = bψ₁ + bψ₂
    bψminus = bψ₂ - bψ₁
    ∂bψplus = ∂bψ₁ + ∂bψ₂
    ∂bψminus = ∂bψ₂ - ∂bψ₁

    Jplus_minus = ψplus_minus * ∂ψplus_minus
    Jplus_plus = ψplus_plus * ∂ψplus_plus
    Jminus_plus = ψminus_plus * ∂ψminus_plus
    Jplus_dagger = ∂bψplus * bψplus
    Jminus_dagger = ∂bψminus * bψminus

    loss =
        (1 // 8) *
        im *
        (
            (Jplus_dagger - Jminus_dagger) * Jplus_minus -
            (Jplus_plus - Jminus_plus) * Jminus_dagger
        )
    return InteractionLagrangian(loss, :γp)
end

function self_energy_diagrams(Σ)
    return [
        diagram for component in (Σ.keldysh, Σ.retarded, Σ.advanced) for
        (diagram, _) in component
    ]
end

function fourier_self_energy_graphs(Σ)
    return [
        graph for component in (Σ.keldysh, Σ.retarded, Σ.advanced) for
        (graph, _) in component
    ]
end

@testset "two-particle irreducibility preserves isolated cut vertices" begin
    c = skeleton_loss_ϕ[Classical]
    q = skeleton_loss_ϕ[Quantum]
    b1 = KC.Bulk(1)
    b2 = KC.Bulk(2)

    sunset = KC.Contraction[(c(b1), bar(q)(b2)), (q(b1), bar(c)(b2)), (c(b1), bar(q)(b2))]
    lollipop = KC.Contraction[(c(b1), bar(q)(b2)), (q(b1), bar(c)(b2)), (c(b2), bar(q)(b2))]

    @test KC.is_irreducible(sunset)
    @test KC.is_two_particle_irreducible(sunset)
    @test KC.is_self_energy_skeleton(sunset)

    @test KC.is_irreducible(lollipop)
    @test !KC.is_two_particle_irreducible(lollipop)
    @test !KC.is_self_energy_skeleton(lollipop)
end

@testset "generated bosonic loss separates 1PI and skeleton self-energies" begin
    L = skeleton_loss_lagrangian()

    G1 = DressedPropagator(L, Val(1), Val(3); simplify=true, _set_reg_to_zero=false)
    ordinary1 = @inferred SelfEnergy(G1)
    skeleton1 = @inferred KC.skeleton_self_energy(G1)
    @test skeleton1.keldysh == ordinary1.keldysh
    @test skeleton1.retarded == ordinary1.retarded
    @test skeleton1.advanced == ordinary1.advanced

    G2 = DressedPropagator(L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false)
    ordinary2 = @inferred SelfEnergy(G2)
    skeleton2 = @inferred KC.skeleton_self_energy(G2)
    ordinary_diagrams = self_energy_diagrams(ordinary2)
    skeleton_diagrams = self_energy_diagrams(skeleton2)

    @test !isempty(skeleton_diagrams)
    @test length(skeleton_diagrams) < length(ordinary_diagrams)
    @test all(diagram -> KC.is_irreducible(KC.contractions(diagram)), ordinary_diagrams)
    @test all(
        diagram -> KC.is_self_energy_skeleton(KC.contractions(diagram)), skeleton_diagrams
    )
    @test any(
        diagram -> !KC.is_two_particle_irreducible(KC.contractions(diagram)),
        ordinary_diagrams,
    )

    ordinary_topologies = Set(only(KC.topology(diagram)) for diagram in ordinary_diagrams)
    skeleton_topologies = Set(only(KC.topology(diagram)) for diagram in skeleton_diagrams)
    @test 2 in ordinary_topologies
    @test 3 in ordinary_topologies
    @test skeleton_topologies == Set([3])
end

@testset "Fourier skeleton selection follows the normal kinetic entry point" begin
    L = skeleton_loss_lagrangian()
    G = DressedPropagator(L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false)
    GF = @inferred fourier_transform(G)
    ordinary = @inferred SelfEnergy(GF)
    skeleton = @inferred KC.skeleton_self_energy(GF)
    ordinary_graphs = fourier_self_energy_graphs(ordinary)
    skeleton_graphs = fourier_self_energy_graphs(skeleton)

    @test !isempty(skeleton_graphs)
    @test length(skeleton_graphs) < length(ordinary_graphs)
    @test all(
        graph -> KC.is_two_particle_irreducible(KC.contractions(graph.coordinate)),
        skeleton_graphs,
    )
    @test any(
        graph -> !KC.is_two_particle_irreducible(KC.contractions(graph.coordinate)),
        ordinary_graphs,
    )

    ordinary_topologies = Set(
        only(KC.topology(graph.coordinate)) for graph in ordinary_graphs
    )
    skeleton_topologies = Set(
        only(KC.topology(graph.coordinate)) for graph in skeleton_graphs
    )
    @test 2 in ordinary_topologies
    @test 3 in ordinary_topologies
    @test skeleton_topologies == Set([3])
end

@testset "fermionic skeleton selector shares the same self-energy machinery" begin
    L = skeleton_pwave_lagrangian()
    G = DressedPropagator(L, Val(1), Val(3); simplify=true, preserve_regularisation=true)
    GF = @inferred fourier_transform(G)
    ordinary = @inferred SelfEnergy(GF)
    skeleton = @inferred KC.skeleton_self_energy(GF)
    @test skeleton == ordinary
end
