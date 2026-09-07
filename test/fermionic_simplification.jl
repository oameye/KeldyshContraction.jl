using KeldyshContraction, Test
import KeldyshContraction as KC

@testset "fermionic default equal-time LO simplification" begin
    @qfields ψsimplify::Fermion χsimplify::Fermion
    ψ₁ = ψsimplify[One]
    χ₂ = χsimplify[Two]

    vertex = @inferred ψ₁ * χ₂ * bar(ψ₁) * bar(χ₂)
    L = @inferred InteractionLagrangian(vertex, :u)
    G = @inferred DressedPropagator(L, Val(1), Val(3); target=ψsimplify)

    @test length(G.retarded) == 1
    @test length(G.keldysh) == 1
    @test isempty(G.advanced)
    @test only(values(G.retarded.diagrams)) == im
    @test only(values(G.keldysh.diagrams)) == im

    retarded_edges = KC.contractions(only(keys(G.retarded.diagrams)))
    @test count(KC.is_retarded, retarded_edges) == 3
    @test count(KC.is_advanced, retarded_edges) == 0
    @test count(KC.is_keldysh, retarded_edges) == 0

    keldysh_edges = KC.contractions(only(keys(G.keldysh.diagrams)))
    @test count(KC.is_retarded, keldysh_edges) == 2
    @test count(KC.is_advanced, keldysh_edges) == 0
    @test count(KC.is_keldysh, keldysh_edges) == 1

    Σ = @inferred SelfEnergy(G)
    @test length(Σ.retarded) == 1
    @test isempty(Σ.keldysh)
    @test isempty(Σ.advanced)
    @test only(values(Σ.retarded.diagrams)) == im

    internal_edge = only(KC.contractions(only(keys(Σ.retarded.diagrams))))
    @test KC.is_retarded(internal_edge)
    @test KC.is_one(internal_edge.out)
    @test KC.is_one(internal_edge.in)
    @test field_family(internal_edge.out) === χsimplify
    @test field_family(internal_edge.in) === χsimplify
end
