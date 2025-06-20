using KeldyshContraction, Test

@qfields c::Destroy(Classical) q::Destroy(Quantum)
elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
L_int = InteractionLagrangian(elasctic2boson)

GF = DressedPropagator(L_int, 2)

Σ = SelfEnergy(GF, 2)

Σk = wigner_transform(Σ)

@testset "reduce to spectral" begin
    using KeldyshContraction: reduce_to_spectral
    tmp = reduce_to_spectral(Σk.keldysh)
    @test length(tmp.diagrams) == 3
    @test Set(abs.(values(tmp.diagrams))) == Set(abs.(values(Σk.keldysh.diagrams)))
end

@testset "reduce to spectral" begin
    using KeldyshContraction: reduce_to_spectral
    tmp = reduce_to_spectral(Σk.keldysh)
    @test length(tmp.diagrams) == 3
    @test Set(abs.(values(tmp.diagrams))) == Set(abs.(values(Σk.keldysh.diagrams)))
end

@testset "keldysh to distribution" begin
    using KeldyshContraction: reduce_to_spectral, kelysh_to_distribution
    tmp = reduce_to_spectral(Σk.keldysh)
    ΣkF = kelysh_to_distribution(tmp)
    iΣkF = im * ΣkF[[3]]
    @test length(tmp.diagrams) == 3
    @test isequal(Set(real.(values(iΣkF.terms))), Set([-1 / 16, -1 / 4, 1 / 2]))
end

@testset "keldysh to distribution" begin
    using KeldyshContraction: reduce_to_spectral, kelysh_to_distribution
    tmp = reduce_to_spectral(Σk.keldysh)
    ΣkF = kelysh_to_distribution(tmp)
    iΣkF = im * ΣkF[[3]]
    @test length(tmp.diagrams) == 3
    @test isequal(Set(real.(values(iΣkF.terms))), Set([-1 / 16, -1 / 4, 1 / 2]))
    number_of_F = Set(map(t -> length(t.momenta), collect(keys(iΣkF.terms))))
    @test isequal(number_of_F, Set([1, 3]))
end

@testset "imaginary part" begin
    using KeldyshContraction: imaginary_part
    imΣr = imaginary_part(Σk.retarded)[[3]]
    @test length(imΣr.terms) == 3
    @test Set(real.(values(imΣr.terms))) == Set([-1 / 16, 1 / 8])

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
    ci = CollisionIntegral(Σk)
    @test isempty(ci.terms[[2]])
    Cint = ci.terms[[3]]
    @test length(Cint) == 6

    @test contains(repr(ci), "Collision integral")
    # @test Set(real.(values(ci.terms))) == Set([-1 / 16, -1 / 4, 1 / 2])
    # number_of_F = Set(map(t -> length(t.momenta), collect(keys(ci.terms))))
    # @test isequal(number_of_F, Set([1, 3]))
end
