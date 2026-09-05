using KeldyshContraction, Test
using KeldyshContraction: Contraction, WickPairing, prepare_args, pairing_sign
import KeldyshContraction as KC

@testset "Wick pairing representation" begin
    @qfields ϕ::Boson
    c, q = ϕ[Classical], ϕ[Quantum]

    args = Field{Boson}[c, q, bar(c), bar(q)]
    destroys, creates = @inferred prepare_args(args, Val(2))
    @test destroys == [c, q]
    @test creates == [bar(q), bar(c)]

    contractions = Contraction{Boson}[Contraction(c, bar(q)), Contraction(q, bar(c))]
    pairing = @inferred WickPairing(contractions, Int8(1), Val(2))
    @test pairing isa WickPairing{Boson,2}
    @test isconcretetype(typeof(pairing))
    @test pairing.sign == Int8(1)
    @test pairing_sign(Boson, (2, 1)) == Int8(1)
end

@testset "field-family compatibility" begin
    @qfields ϕ::Boson χ::Boson
    c, q = ϕ[Classical], ϕ[Quantum]
    χq = χ[Quantum]

    @test KC.contraction_compatible(Boson, c, bar(q))
    @test !KC.contraction_compatible(Boson, c, bar(χq))
end
