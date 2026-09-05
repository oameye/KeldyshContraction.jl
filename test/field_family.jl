using KeldyshContraction, Test
using KeldyshContraction: FieldFamily, FieldIndex, FieldIndices, field_family, In, Out
using KeldyshContraction: Regularisation.Plus as Plus

@testset "FieldFamily identity" begin
    indices = FieldIndices(FieldIndex(:spin, 1), FieldIndex(:species, 2))
    family = FieldFamily{Boson}(:ψ, indices)

    c = @inferred Field(family, Classical)
    q = @inferred Field(family, Quantum)

    @test typeof(c) === Field{Boson}
    @test typeof(q) === Field{Boson}
    @test isconcretetype(typeof(family))
    @test all(isconcretetype, fieldtypes(typeof(family)))

    @test @inferred(field_family(c)) == family
    @test @inferred(field_family(q)) == family
    @test field_family(bar(q)(In())(Plus)) == family
    @test field_family(c(Out())) == family
    @test hash(field_family(c)) == hash(field_family(q))

    @test c != q
    @test FieldFamily{Boson}(:χ, indices) != family
    @test FieldFamily{Boson}(:ψ, FieldIndices(FieldIndex(:spin, 2))) != family
end
