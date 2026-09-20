using KeldyshContraction, Test
import KeldyshContraction as KC

@testset "Direct self-energy generation" begin
    @qfields direct_Σ_ϕ::Boson
    c, q = direct_Σ_ϕ[Classical], direct_Σ_ϕ[Quantum]
    elastic = -(0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2))
    L = InteractionLagrangian(elastic)

    for (order, edges) in ((Val(2), Val(5)), (Val(3), Val(7)))
        reference = SelfEnergy(DressedPropagator(L, order, edges))
        direct = KC._direct_self_energy(L, order, edges)

        @test direct.keldysh == reference.keldysh
        @test direct.retarded == reference.retarded
        @test direct.advanced == reference.advanced
        @test KC.parameters(direct) == KC.parameters(reference)
        @test KC.target_family(direct) == KC.target_family(reference)
        @test KC.order(direct) == KC.order(reference)
    end
end
