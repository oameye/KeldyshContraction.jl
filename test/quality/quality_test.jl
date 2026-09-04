using KeldyshContraction
using Test
using Aqua
using CheckConcreteStructs: all_concrete
using ExplicitImports

@testset "Concretely typed" begin
    all_concrete(KeldyshContraction.Field{KeldyshContraction.Boson})
    all_concrete(KeldyshContraction.FieldIndex)
    all_concrete(KeldyshContraction.FieldIndices)
    all_concrete(KeldyshContraction.QMul)
    all_concrete(KeldyshContraction.QAdd)
    all_concrete(KeldyshContraction.InteractionLagrangian)

    all_concrete(KeldyshContraction.Edge)
    all_concrete(KeldyshContraction.Diagram{5,3})
    all_concrete(KeldyshContraction.Diagrams{5,3})
    all_concrete(KeldyshContraction.DressedPropagator)
    all_concrete(KeldyshContraction.SelfEnergy)
end

@testset "ExplicitImports" begin
    allow_unanalyzable = (
        KeldyshContraction.Regularisation,
        KeldyshContraction.PropagatorType,
        KeldyshContraction.KeldyshIndex,
        KeldyshContraction.Orientation,
    )

    @test check_no_implicit_imports(KeldyshContraction; allow_unanalyzable) == nothing
    @test check_all_explicit_imports_via_owners(KeldyshContraction) == nothing
    @test check_all_explicit_imports_are_public(KeldyshContraction) == nothing
    @test check_no_stale_explicit_imports(KeldyshContraction; allow_unanalyzable) == nothing
    @test check_all_qualified_accesses_via_owners(KeldyshContraction) == nothing
    @test check_no_self_qualified_accesses(KeldyshContraction) == nothing
end

@testset "best practices" begin
    Aqua.test_ambiguities([KeldyshContraction]; broken=false)
    Aqua.test_all(KeldyshContraction; ambiguities=false)
end
