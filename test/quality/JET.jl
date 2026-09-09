using KeldyshContraction
using KeldyshContraction:
    Bulk,
    Contraction,
    Diagram,
    FourierDiagram,
    In,
    MomentumBasis,
    Out,
    RoutingEdge,
    basis_momentum,
    edge_momenta,
    kinematic_factor,
    lower_fourier_derivatives,
    matrix,
    momentum_basis,
    momentum_routing,
    routing_matrix
using Test
using JET

@qfields jet_ψ_raw::Fermion jet_χ_raw::Fermion
@qfields jet_ϕ_raw::Boson
const jet_ψ = jet_ψ_raw
const jet_χ = jet_χ_raw
const jet_ψ₁ = jet_ψ[One]
const jet_χ₂ = jet_χ[Two]
const jet_ϕc = jet_ϕ_raw[Classical]
const jet_ϕq = jet_ϕ_raw[Quantum]

function fermionic_jet_workload()
    interaction = jet_ψ₁ * jet_χ₂ * bar(jet_ψ₁) * bar(jet_χ₂)
    L = InteractionLagrangian(interaction)
    G = DressedPropagator(L, Val(1), Val(3); target=jet_ψ)
    Σ = SelfEnergy(G)
    return matrix(G), matrix(Σ)
end

function fermionic_lagrangian_sum_jet_workload()
    interaction = jet_ψ₁ * jet_χ₂ * bar(jet_ψ₁) * bar(jet_χ₂)
    L_u = InteractionLagrangian(interaction, :u)
    L_v = InteractionLagrangian(interaction, :v)
    G = DressedPropagator(L_u + L_v, Val(1), Val(3); target=jet_ψ)
    Σ = SelfEnergy(G)
    return matrix(G[:u]), matrix(Σ[:v])
end

function derivative_fermionic_jet_workload()
    ∂xχ₂ = partial(jet_χ₂, :x)
    interaction = jet_ψ₁ * ∂xχ₂ * bar(jet_ψ₁) * bar(∂xχ₂)
    L = InteractionLagrangian(interaction, :d)
    G = DressedPropagator(L, Val(1), Val(3); target=jet_ψ, simplify=false)
    Σ = SelfEnergy(G)
    return derivatives(∂xχ₂), matrix(G), matrix(Σ)
end

function momentum_routing_jet_workload()
    basis = MomentumBasis(2)
    linear = basis_momentum(basis, 1) - 2 * basis_momentum(basis, 2)
    edges = RoutingEdge[
        RoutingEdge(1, 2), RoutingEdge(2, 3), RoutingEdge(3, 1), RoutingEdge(1, 3)
    ]
    routing = momentum_routing(edges)
    return linear, length(routing.basis), routing_matrix(routing)
end

function fourier_diagram_jet_workload()
    contractions = Contraction{Boson}[
        Contraction(jet_ϕc(Out()), bar(jet_ϕq)(Bulk(1))),
        Contraction(jet_ϕc(Bulk(1)), bar(jet_ϕq)(Bulk(1))),
        Contraction(jet_ϕc(Bulk(1)), bar(jet_ϕq)(In())),
    ]
    diagram = Diagram(contractions, Val(3), Val(0))
    routed = FourierDiagram(diagram)
    return length(momentum_basis(routed)), edge_momenta(routed)
end

function fourier_derivative_jet_workload()
    ∂xϕc = partial(jet_ϕc, :x)
    contractions = Contraction{Boson}[
        Contraction(∂xϕc(Out()), bar(jet_ϕq)(Bulk(1))),
        Contraction(jet_ϕc(Bulk(1)), bar(jet_ϕq)(Bulk(1))),
        Contraction(jet_ϕc(Bulk(1)), bar(jet_ϕq)(In())),
    ]
    diagram = Diagram(contractions, Val(3), Val(0))
    routed = FourierDiagram(diagram)
    lowered = lower_fourier_derivatives(routed)
    return kinematic_factor(lowered)
end

function fourier_pwave_jet_workload()
    ψ₂ = jet_ψ[Two]
    ∂xψ₂ = partial(ψ₂, :x)
    interaction = jet_ψ₁ * ∂xψ₂ * bar(jet_ψ₁) * bar(∂xψ₂)
    L = InteractionLagrangian(interaction, :d)
    G = DressedPropagator(L, Val(1), Val(3); target=jet_ψ, simplify=false)
    Gk = fourier_transform(G)
    Σk = SelfEnergy(Gk)
    return matrix(Gk), matrix(Σk)
end

@static if isempty(VERSION.prerelease)
    @testset "JET report_package" begin
        rep = JET.report_package(KeldyshContraction; target_modules=(KeldyshContraction,))
        @show rep
        @test isempty(JET.get_reports(rep))
    end

    @testset "JET fermionic public workload" begin
        JET.@test_opt target_modules=(KeldyshContraction,) fermionic_jet_workload()
    end

    @testset "JET fermionic LagrangianSum workload" begin
        JET.@test_opt target_modules=(KeldyshContraction,) fermionic_lagrangian_sum_jet_workload()
    end

    @testset "JET derivative fermionic workload" begin
        JET.@test_opt target_modules=(KeldyshContraction,) derivative_fermionic_jet_workload()
    end

    @testset "JET exact momentum-routing workload" begin
        JET.@test_opt target_modules=(KeldyshContraction,) momentum_routing_jet_workload()
    end

    @testset "JET Fourier-diagram routing workload" begin
        JET.@test_opt target_modules=(KeldyshContraction,) fourier_diagram_jet_workload()
    end

    @testset "JET Fourier derivative-lowering workload" begin
        JET.@test_opt target_modules=(KeldyshContraction,) fourier_derivative_jet_workload()
    end

    @testset "JET p-wave Fourier self-energy workload" begin
        JET.@test_opt target_modules=(KeldyshContraction,) fourier_pwave_jet_workload()
    end
end
