using KeldyshContraction, Test

@testset "Concretely typed" begin
    using KeldyshContraction
    using CheckConcreteStructs

    all_concrete(KeldyshContraction.QMul)
    all_concrete(KeldyshContraction.QAdd)
    all_concrete(KeldyshContraction.Destroy)
    all_concrete(KeldyshContraction.Create)
    all_concrete(KeldyshContraction.Propagator)
end

@testset "Code linting" begin
    using JET
    JET.test_package(KeldyshContraction; target_defined_modules=true)
end

@testset "ExplicitImports" begin
    using ExplicitImports
    @test check_no_stale_explicit_imports(KeldyshContraction) == nothing
    @test check_all_explicit_imports_via_owners(KeldyshContraction) == nothing
end

@testset "best practices" begin
    using Aqua

    Aqua.test_ambiguities([KeldyshContraction], broken=false)
    Aqua.test_all(KeldyshContraction; ambiguities=false)
end

@testset "propagator" begin
    include("propagator.jl")
end

@testset "two_boson_loss" begin
    include("two_boson_loss.jl")
end
