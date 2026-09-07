using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields ϕ::Boson
c, q = ϕ[Classical], ϕ[Quantum]
elasctic2boson = -(
    0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2)
)
L_int = InteractionLagrangian(elasctic2boson)

GF = @inferred DressedPropagator(L_int, Val(2), Val(5))
Σ = @inferred SelfEnergy(GF)
Σk = @inferred wigner_transform(Σ)

@testset "distribution coefficient representation" begin
    using KeldyshContraction: BosonicDistributionTerm, BosonicDistributions, Momenta

    @test isempty(@inferred(BosonicDistributions()))

    term = BosonicDistributionTerm([Momenta(0)])
    real_distribution = BosonicDistributions{Float64}(Dict(term => 1.0))
    complex_distribution = @inferred im * real_distribution

    @test complex_distribution isa BosonicDistributions{ComplexF64}
    @test valtype(typeof(complex_distribution.terms)) === ComplexF64
end

@testset "reduce to spectral" begin
    using KeldyshContraction: reduce_to_spectral
    tmp = @inferred reduce_to_spectral(Σk.keldysh)
    @test typeof(tmp) === typeof(Σk.keldysh)
    @test length(tmp.diagrams) == 3
    @test Set(abs.(values(tmp.diagrams))) == Set(abs.(values(Σk.keldysh.diagrams)))
end

@testset "reduce to spectral without simplification" begin
    using KeldyshContraction: reduce_to_spectral, Bulk, Diagram, Diagrams, Edge
    d = Diagram(
        [Edge(c(Bulk(1)), bar(c)(Bulk(2))), Edge(c(Bulk(3)), bar(c)(Bulk(4)))],
        Val(2),
        Val(6),
    )
    ds = Diagrams([d], Complex{Rational{Int}}(1.0))
    @test isequal(@inferred(reduce_to_spectral(ds)), ds)
end

@testset "reduce to spectral duplicate check" begin
    using KeldyshContraction: reduce_to_spectral
    tmp = @inferred reduce_to_spectral(Σk.keldysh)
    @test length(tmp.diagrams) == 3
    @test Set(abs.(values(tmp.diagrams))) == Set(abs.(values(Σk.keldysh.diagrams)))
end

@testset "keldysh to distribution" begin
    using KeldyshContraction:
        BosonicDistributionTerm,
        BosonicDistributions,
        reduce_to_spectral,
        kelysh_to_distribution
    tmp = @inferred reduce_to_spectral(Σk.keldysh)
    ΣkF = @inferred kelysh_to_distribution(tmp)
    @test ΣkF isa Dict{KC.FixedVector{1,Int},BosonicDistributions{ComplexF64}}
    iΣkF = @inferred im * ΣkF[[3]]
    @test iΣkF isa BosonicDistributions{ComplexF64}
    @test length(tmp.diagrams) == 3
    @test isequal(Set(real.(values(iΣkF.terms))), Set([-1 / 16, -1 / 4, 1 / 2]))

    diagram = first(keys(tmp.diagrams))
    single_distribution = @inferred kelysh_to_distribution(diagram)
    @test first(single_distribution) isa BosonicDistributionTerm
    @test last(single_distribution) isa Number
end

@testset "keldysh to distribution multiplicity" begin
    using KeldyshContraction: reduce_to_spectral, kelysh_to_distribution
    tmp = @inferred reduce_to_spectral(Σk.keldysh)
    ΣkF = @inferred kelysh_to_distribution(tmp)
    iΣkF = @inferred im * ΣkF[[3]]
    @test length(tmp.diagrams) == 3
    @test isequal(Set(real.(values(iΣkF.terms))), Set([-1 / 16, -1 / 4, 1 / 2]))
    number_of_F = Set(map(t -> length(t.momenta), collect(keys(iΣkF.terms))))
    @test isequal(number_of_F, Set([1, 3]))
end

@testset "imaginary part" begin
    using KeldyshContraction: BosonicDistributionTerm, BosonicDistributions, imaginary_part
    imΣr_all = @inferred imaginary_part(Σk.retarded)
    @test imΣr_all isa Dict{KC.FixedVector{1,Int},BosonicDistributions{ComplexF64}}
    imΣr = imΣr_all[[3]]
    @test length(imΣr.terms) == 3
    @test Set(real.(values(imΣr.terms))) == Set([-1 / 16, 1 / 8])

    diagram = first(keys(Σk.retarded.diagrams))
    single_imaginary = @inferred imaginary_part(diagram)
    @test first(single_imaginary) isa BosonicDistributionTerm
    @test last(single_imaginary) isa Number

    imΣr_string = repr(imΣr)
    @test contains(imΣr_string, "-0.0625*F(q₁)*F(q₂)")
    @test contains(imΣr_string, "0.125*F(q₁)*F(q₁ + q₂ - k)")
    @test contains(imΣr_string, "-0.0625")

    @test isequal(
        Set(repr.(collect(keys(imΣr.terms)))),
        Set(["F(q₁)*F(q₂)", "", "F(q₁)*F(q₁ + q₂ - k)"]),
    )
end

@testset "Collision Integral" begin
    using KeldyshContraction: CollisionIntegral
    ci = @inferred CollisionIntegral(Σk)
    @test ci isa CollisionIntegral{ComplexF64,1}
    @test isempty(ci.terms[[2]])
    Cint = ci.terms[[3]]
    @test length(Cint) == 6

    @test contains(repr(ci), "Collision integral")
end
