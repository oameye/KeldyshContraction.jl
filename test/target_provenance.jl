using KeldyshContraction, Test

@qfields provenance_ψ::Fermion provenance_χ::Fermion

@testset "external target family survives the two-point pipeline" begin
    ψ₁ = provenance_ψ[One]
    χ₂ = provenance_χ[Two]
    ∂xχ₂ = partial(χ₂, :x)
    L = InteractionLagrangian((1 // 2) * ψ₁ * ∂xχ₂ * bar(ψ₁) * bar(∂xχ₂), :d)

    G = @inferred DressedPropagator(L, Val(1), Val(3); target=provenance_ψ, simplify=false)
    @test @inferred(target_family(G)) === provenance_ψ

    GF = @inferred fourier_transform(G)
    @test @inferred(target_family(GF)) === provenance_ψ

    ΣF = @inferred SelfEnergy(GF)
    @test @inferred(target_family(ΣF)) === provenance_ψ

    ΣW = @inferred wigner_transform(ΣF; gradient_order=Val(0))
    @test @inferred(target_family(ΣW)) === provenance_ψ

    KΣ = @inferred kinetic_expression(ΣW)
    @test @inferred(target_family(KΣ)) === provenance_ψ
end
