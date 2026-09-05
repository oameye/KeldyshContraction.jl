using KeldyshContraction, Test
using KeldyshContraction: FieldIndex, FieldIndices, In, Out
using KeldyshContraction: Regularisation.Plus as Plus

@testset "FieldFamily identity" begin
    indices = FieldIndices(FieldIndex(:spin, 1), FieldIndex(:species, 2))
    family = FieldFamily{Boson}(:ψ, indices)

    c = @inferred family[Classical]
    q = @inferred family[Quantum]

    @test typeof(c) === Field{Boson}
    @test typeof(q) === Field{Boson}
    @test isconcretetype(typeof(family))
    @test all(isconcretetype, fieldtypes(typeof(family)))

    @test @inferred(field_family(c)) === family
    @test @inferred(field_family(q)) === family
    @test field_family(bar(q)(In())(Plus)) === family
    @test field_family(c(Out())) === family
    @test hash(field_family(c)) == hash(field_family(q))

    @test c != q
    @test FieldFamily{Boson}(:χ, indices) != family
    @test FieldFamily{Boson}(:ψ, FieldIndices(FieldIndex(:spin, 2))) != family
end

@testset "field-family declaration" begin
    @qfields ϕ::Boson χ::Boson

    c = @inferred ϕ[Classical]
    q = @inferred ϕ[Quantum]
    χc = @inferred χ[Classical]

    @test ϕ isa FieldFamily{Boson}
    @test χ isa FieldFamily{Boson}
    @test field_family(c) === ϕ
    @test field_family(q) === ϕ
    @test field_family(χc) === χ
    @test field_family(c) != field_family(χc)
end

@testset "interaction field families" begin
    @qfields ϕ::Boson χ::Boson
    c, q = ϕ[Classical], ϕ[Quantum]
    χc, χq = χ[Classical], χ[Quantum]

    Lϕ = InteractionLagrangian(bar(c) * q + bar(q) * c)
    @test @inferred(field_families(Lϕ)) == [ϕ]

    Lϕχ = InteractionLagrangian(
        bar(c) * q + bar(q) * c + bar(χc) * χq + bar(χq) * χc
    )
    @test Set(@inferred(field_families(Lϕχ))) == Set([ϕ, χ])
end
