using KeldyshContraction, Test

@testset "Concretely typed" begin
    using KeldyshContraction
    using CheckConcreteStructs

    all_concrete(KeldyshContraction.QMul)
    all_concrete(KeldyshContraction.QAdd)
    all_concrete(KeldyshContraction.Destroy)
    all_concrete(KeldyshContraction.Create)
    all_concrete(KeldyshContraction.InteractionLagrangian)

    all_concrete(KeldyshContraction.Edge)
    all_concrete(KeldyshContraction.Diagram)
    all_concrete(KeldyshContraction.Diagrams)
    all_concrete(KeldyshContraction.DressedPropagator)
    all_concrete(KeldyshContraction.SelfEnergy)
end

if isempty(VERSION.prerelease)
    @testset "Code linting" begin
        using JET
        JET.test_package(KeldyshContraction; target_defined_modules=true)
        rep = report_package("KeldyshContraction")
        @show rep
        @test length(JET.get_reports(rep)) <= 24
        @test_broken length(JET.get_reports(rep)) == 0
    end
end

@testset "ExplicitImports" begin
    using ExplicitImports
    @test check_no_stale_explicit_imports(KeldyshContraction) == nothing
    @test check_all_explicit_imports_via_owners(KeldyshContraction) == nothing
end

@testset "best practices" begin
    using Aqua

    Aqua.test_ambiguities([KeldyshContraction]; broken=false)
    Aqua.test_all(KeldyshContraction; ambiguities=false)
end

@testset "show methods" begin
    include("show_methods.jl")
end

@testset "QField" begin
    include("QField.jl")
end

@testset "propagator" begin
    include("propagator.jl")
end

@testset "tadpole reguralisation" begin
    include("reguralise.jl")
end

@testset "diagram" begin
    include("diagram.jl")
end

@testset "self-energy" begin
    include("self_energy.jl")
end

@testset "two boson loss" begin
    include("two_boson_loss.jl")
end

@testset "elastic two body" begin
    include("elastic_two_body.jl")
end

@testset "self interaction" begin
    include("self_interaction.jl")
end

@testset "second order" begin
    include("second_order.jl")
end
