using KeldyshContraction, Test
import KeldyshContraction as KC
using KeldyshContraction:
    Bulk, Contraction, Diagram, In, LinearMomentum, MomentumVariable, Out

function wigner_recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        wigner_recursively_concrete(eltype(T), seen) || return false
    elseif T <: AbstractDict
        wigner_recursively_concrete(keytype(T), seen) || return false
        wigner_recursively_concrete(valtype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        wigner_recursively_concrete(FT, seen) || return false
    end
    return true
end

function assert_wigner_preserves_fourier(wigner, fourier)
    @test length(wigner) == length(fourier)
    @test wigner_context(wigner) == HomogeneousWignerContext()
    @test gradient_order(wigner) == Val(0)

    fourier_by_graph = Dict(graph => contributions for (graph, contributions) in fourier)
    for (wigner_graph, wigner_contributions) in wigner
        fourier_graph = wigner_graph.fourier
        @test haskey(fourier_by_graph, fourier_graph)
        @test gradient_order(wigner_graph) == Val(0)
        @test KC.external_momentum_count(wigner_graph) == 1
        @test external_wigner_momentum(wigner_graph) == KC.momentum_basis(fourier_graph)[1]
        @test KC.momentum_basis(wigner_graph) == KC.momentum_basis(fourier_graph)
        @test KC.edge_momenta(wigner_graph) == KC.edge_momenta(fourier_graph)
        @test KC.loop_momentum_count(wigner_graph) == KC.loop_momentum_count(fourier_graph)
        @test KC.coordinate_diagram(wigner_graph) == KC.coordinate_diagram(fourier_graph)

        fourier_contributions = fourier_by_graph[fourier_graph]
        @test length(wigner_contributions) == length(fourier_contributions)
        for (wc, fc) in zip(wigner_contributions, fourier_contributions)
            @test wc.coefficient == fc.coefficient
            @test wc.kinematic == fc.kinematic
        end
    end
end

@qfields wigner_pwave_ψ::Fermion

@testset "homogeneous Wigner self-energy preserves exact Fourier data" begin
    ψ₁ = wigner_pwave_ψ[One]
    ψ₂ = wigner_pwave_ψ[Two]
    ∂xψ₂ = partial(ψ₂, :x)
    vertex = @inferred ψ₁ * ∂xψ₂ * bar(ψ₁) * bar(∂xψ₂)

    L = @inferred InteractionLagrangian(vertex, :γ)
    G = @inferred DressedPropagator(L, Val(1), Val(3); simplify=false)
    Gk = @inferred fourier_transform(G)
    Σk = @inferred SelfEnergy(Gk)
    ΣW = @inferred wigner_transform(Σk; gradient_order=Val(0))
    ΣW_default = @inferred wigner_transform(Σk)

    @test ΣW isa WignerSelfEnergy
    @test ΣW_default == ΣW
    @test KC.statistics(ΣW) === Fermion
    @test KC.order(ΣW) == 1
    @test parameters(ΣW) == parameters(Σk)
    @test gradient_order(ΣW) == Val(0)
    @test wigner_recursively_concrete(typeof(ΣW))
    @test wigner_context(ΣW) == HomogeneousWignerContext()
    @test wigner_context(ΣW).center_coordinate === nothing

    assert_wigner_preserves_fourier(ΣW.retarded, Σk.retarded)
    assert_wigner_preserves_fourier(ΣW.keldysh, Σk.keldysh)
    assert_wigner_preserves_fourier(ΣW.advanced, Σk.advanced)

    ΣW_again = @inferred wigner_transform(Σk; gradient_order=Val(0))
    @test ΣW == ΣW_again
    @test isequal(ΣW, ΣW_again)
    @test hash(ΣW) == hash(ΣW_again)

    @test_throws ArgumentError wigner_transform(Σk; gradient_order=Val(1))
    @test_throws ArgumentError wigner_transform(@inferred(SelfEnergy(G)))
end

@testset "homogeneous Wigner dressed propagator public path" begin
    ψ₁ = wigner_pwave_ψ[One]
    ψ₂ = wigner_pwave_ψ[Two]
    ∂xψ₂ = partial(ψ₂, :x)
    L = @inferred InteractionLagrangian(ψ₁ * ∂xψ₂ * bar(ψ₁) * bar(∂xψ₂), :γ)
    G = @inferred DressedPropagator(L, Val(1), Val(3); simplify=false)
    Gk = @inferred fourier_transform(G)

    GW = @inferred wigner_transform(Gk; gradient_order=Val(0))
    GW_default = @inferred wigner_transform(Gk)
    GW_coordinate = @inferred wigner_transform(G; gradient_order=Val(0))

    @test GW isa WignerDressedPropagator
    @test GW_default == GW
    @test GW == GW_coordinate
    @test KC.statistics(GW) === Fermion
    @test KC.order(GW) == 1
    @test parameters(GW) == parameters(Gk)
    @test gradient_order(GW) == Val(0)
    @test wigner_recursively_concrete(typeof(GW))

    assert_wigner_preserves_fourier(GW.retarded, Gk.retarded)
    assert_wigner_preserves_fourier(GW.keldysh, Gk.keldysh)
    assert_wigner_preserves_fourier(GW.advanced, Gk.advanced)

    @test KC.matrix(GW)[1, 1] == GW.retarded
    @test KC.matrix(GW)[1, 2] == GW.keldysh
    @test KC.matrix(GW)[2, 2] == GW.advanced
    @test_throws ArgumentError wigner_transform(Gk; gradient_order=Val(2))
end

@qfields wigner_boson_ϕ::Boson
@qfields wigner_fermion_ψ::Fermion

function routing_parity_diagram(::Type{Boson})
    c = wigner_boson_ϕ[Classical]
    q = wigner_boson_ϕ[Quantum]
    contractions = Contraction{Boson}[
        Contraction(c(Out()), bar(q)(Bulk(1))),
        Contraction(c(Bulk(1)), bar(q)(Bulk(1))),
        Contraction(c(Bulk(1)), bar(q)(In())),
    ]
    return Diagram(contractions, Val(3), Val(0))
end

function routing_parity_diagram(::Type{Fermion})
    one = wigner_fermion_ψ[One]
    two = wigner_fermion_ψ[Two]
    contractions = Contraction{Fermion}[
        Contraction(one(Out()), bar(two)(Bulk(1))),
        Contraction(one(Bulk(1)), bar(two)(Bulk(1))),
        Contraction(one(Bulk(1)), bar(two)(In())),
    ]
    return Diagram(contractions, Val(3), Val(0))
end

@testset "Wigner routing is statistics-neutral" begin
    boson_fourier = @inferred fourier_transform(routing_parity_diagram(Boson))
    fermion_fourier = @inferred fourier_transform(routing_parity_diagram(Fermion))
    boson_wigner = @inferred wigner_transform(boson_fourier; gradient_order=Val(0))
    fermion_wigner = @inferred wigner_transform(fermion_fourier; gradient_order=Val(0))

    @test KC.momentum_basis(boson_wigner) == KC.momentum_basis(fermion_wigner)
    @test KC.edge_momenta(boson_wigner) == KC.edge_momenta(fermion_wigner)
    @test KC.external_momentum_count(boson_wigner) == 1
    @test KC.external_momentum_count(fermion_wigner) == 1
    @test KC.loop_momentum_count(boson_wigner) == KC.loop_momentum_count(fermion_wigner)
    @test external_wigner_momentum(boson_wigner) == MomentumVariable(1)
    @test external_wigner_momentum(fermion_wigner) == MomentumVariable(1)
    @test KC.statistics(boson_wigner) === Boson
    @test KC.statistics(fermion_wigner) === Fermion
end

@qfields wigner_order2_ϕ::Boson

@testset "representative second-order homogeneous Wigner self-energy" begin
    c = wigner_order2_ϕ[Classical]
    q = wigner_order2_ϕ[Quantum]
    interaction = -(
        1 // 2 * (c^2 + q^2) * bar(c) * bar(q) + 1 // 2 * c * q * (bar(c)^2 + bar(q)^2)
    )
    L = InteractionLagrangian(interaction)
    G = DressedPropagator(L, Val(2), Val(5); simplify=false)
    Σk = SelfEnergy(fourier_transform(G))
    ΣW = @inferred wigner_transform(Σk; gradient_order=Val(0))

    @test ΣW isa WignerSelfEnergy
    @test KC.statistics(ΣW) === Boson
    @test KC.order(ΣW) == 2
    @test wigner_recursively_concrete(typeof(ΣW))
    assert_wigner_preserves_fourier(ΣW.retarded, Σk.retarded)
    assert_wigner_preserves_fourier(ΣW.keldysh, Σk.keldysh)
    assert_wigner_preserves_fourier(ΣW.advanced, Σk.advanced)

    M = KC.matrix(ΣW)
    @test iszero(M[1, 1])
    @test M[1, 2] == ΣW.advanced
    @test M[2, 1] == ΣW.retarded
    @test M[2, 2] == ΣW.keldysh
end

function higher_order_wigner_diagram()
    c = wigner_boson_ϕ[Classical]
    q = wigner_boson_ϕ[Quantum]
    contractions = Contraction{Boson}[
        Contraction(c(Out()), bar(q)(Bulk(1))),
        Contraction(c(Bulk(2)), bar(q)(Bulk(1))),
        Contraction(c(Bulk(2)), bar(q)(Bulk(1))),
        Contraction(c(Bulk(3)), bar(q)(Bulk(2))),
        Contraction(c(Bulk(3)), bar(q)(Bulk(2))),
        Contraction(c(Bulk(1)), bar(q)(Bulk(3))),
        Contraction(c(Bulk(3)), bar(q)(In())),
    ]
    return Diagram(contractions, Val(7), Val(3))
end

@testset "higher-order loop routing survives Wigner representation" begin
    fourier = @inferred fourier_transform(higher_order_wigner_diagram())
    wigner = @inferred wigner_transform(fourier; gradient_order=Val(0))

    @test KC.loop_momentum_count(fourier) == 3
    @test KC.loop_momentum_count(wigner) == 3
    @test KC.fourier_diagram(wigner) == fourier
    @test KC.momentum_basis(wigner) == KC.momentum_basis(fourier)
    @test KC.edge_momenta(wigner) == KC.edge_momenta(fourier)
    @test external_wigner_momentum(wigner) == KC.momentum_basis(fourier)[1]
    @test gradient_order(wigner) == Val(0)
    @test wigner_recursively_concrete(typeof(wigner))
end
