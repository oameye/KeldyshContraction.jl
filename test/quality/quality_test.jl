using KeldyshContraction
using Test
using Aqua
using CheckConcreteStructs: all_concrete
using ExplicitImports

@testset "Concretely typed" begin
    all_concrete(KeldyshContraction.FieldFamily{KeldyshContraction.Boson})
    all_concrete(KeldyshContraction.Field{KeldyshContraction.Boson})
    all_concrete(KeldyshContraction.FieldIndex)
    all_concrete(KeldyshContraction.FieldIndices)
    all_concrete(KeldyshContraction.ParameterPower)
    all_concrete(KeldyshContraction.ParameterMonomial)
    all_concrete(KeldyshContraction.QMul)
    all_concrete(KeldyshContraction.QAdd)
    all_concrete(KeldyshContraction.InteractionLagrangian)

    all_concrete(KeldyshContraction.Momentum)
    all_concrete(KeldyshContraction.Momenta)
    all_concrete(KeldyshContraction.Contraction{KeldyshContraction.Boson})
    all_concrete(KeldyshContraction.Edge{KeldyshContraction.Boson})
    all_concrete(KeldyshContraction.Diagram{KeldyshContraction.Boson,5,3})
    all_concrete(
        KeldyshContraction.Diagrams{
            KeldyshContraction.ComplexRationals,KeldyshContraction.Boson,5,3
        },
    )
    all_concrete(
        KeldyshContraction.DressedPropagator{
            KeldyshContraction.ComplexRationals,KeldyshContraction.Boson,2,5,1
        },
    )
    all_concrete(
        KeldyshContraction.SelfEnergy{
            KeldyshContraction.ComplexRationals,KeldyshContraction.Boson,2,3,1
        },
    )

    all_concrete(KeldyshContraction.BosonicDistributionTerm)
    all_concrete(KeldyshContraction.BosonicDistributions{ComplexF64})
    all_concrete(KeldyshContraction.CollisionIntegral{ComplexF64,1})
end

@testset "No DispatchDoctor exemptions" begin
    src = joinpath(pkgdir(KeldyshContraction), "src")
    for (root, _, files) in walkdir(src)
        for file in files
            endswith(file, ".jl") || continue
            path = joinpath(root, file)
            @test !occursin("@unstable", read(path, String))
        end
    end
end

@testset "ExplicitImports" begin
    allow_unanalyzable = (
        KeldyshContraction.Regularisation,
        KeldyshContraction.PropagatorType,
        KeldyshContraction.KeldyshIndex,
        KeldyshContraction.Orientation,
        KeldyshContraction.IndexKind,
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
