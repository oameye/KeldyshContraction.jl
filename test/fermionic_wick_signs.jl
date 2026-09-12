using KeldyshContraction, Test
import KeldyshContraction as KC

@testset "fermionic four-field Wick signs" begin
    @qfields ψpair::Fermion
    ψ₁ = ψpair[One](KC.Bulk(1))
    ψ₂ = ψpair[Two](KC.Bulk(1))
    bψ₂a = bar(ψpair[Two])(KC.Bulk(2))
    bψ₂b = bar(ψpair[Two])(KC.Bulk(3))

    term = @inferred ψ₁ * ψ₂ * bψ₂a * bψ₂b
    args = copy(KC.fields(term))
    pairings = KC._wick_contraction(args, Val(2); regularise=false)

    @test length(pairings) == 2
    @test all(p -> p isa KC.WickPairing{Fermion,2}, pairings)
    @test sort(Int[p.sign for p in pairings]) == [-1, 1]

    merged = KC._wick_contraction(term, Val(2); regularise=false, simplify=false)
    @test iszero(merged)
end
