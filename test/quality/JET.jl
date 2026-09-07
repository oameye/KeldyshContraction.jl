using KeldyshContraction
using KeldyshContraction: matrix
using Test
using JET

@qfields jet_ψ_raw::Fermion jet_χ_raw::Fermion
const jet_ψ = jet_ψ_raw
const jet_χ = jet_χ_raw
const jet_ψ₁ = jet_ψ[One]
const jet_χ₂ = jet_χ[Two]

function fermionic_jet_workload()
    interaction = jet_ψ₁ * jet_χ₂ * bar(jet_ψ₁) * bar(jet_χ₂)
    L = InteractionLagrangian(interaction)
    G = DressedPropagator(L, Val(1), Val(3); target=jet_ψ)
    Σ = SelfEnergy(G)
    return matrix(G), matrix(Σ)
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
end
