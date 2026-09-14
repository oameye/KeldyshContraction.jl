using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields collision_reduction_ϕ::Boson
@qfields collision_reduction_ψ::Fermion

@testset "statistical to occupation algebra" begin
    k = KC.LinearMomentum([1, 0])
    q = KC.LinearMomentum([0, 1])
    Fk = KC.StatisticalAtom(collision_reduction_ϕ, k)
    Fq = KC.StatisticalAtom(collision_reduction_ϕ, q)
    one = KC.StatisticalMonomial{Boson}()
    mk = KC.StatisticalMonomial([Fk])
    mq = KC.StatisticalMonomial([Fq])
    mkq = KC.StatisticalMonomial([Fk, Fq])
    C = KC.ComplexRationals
    polynomial = KC.StatisticalPolynomial{C,Boson}([
        one => C(-2), mq => C(2), mk => C(2), mkq => C(-2)
    ])
    occupation = @inferred KC.occupation_collision_polynomial(polynomial)
    @test length(occupation) == 1
    monomial, coefficient = only(occupation)
    @test length(monomial) == 2
    @test coefficient == -4 // 1

    kf = KC.LinearMomentum([1])
    Ff = KC.StatisticalAtom(collision_reduction_ψ, kf)
    fermion = KC.StatisticalPolynomial(Ff, KC.ComplexRationals(1))
    nf = @inferred KC.occupation_collision_polynomial(fermion)
    @test length(nf) == 2
    @test KC.occupation_statistics_sign(Fermion) == -1
    @test KC.occupation_collision_factor(Fermion) == -1 // 2
end

function generated_collision_reduction_fixture()
    c = collision_reduction_ϕ[Classical]
    q = collision_reduction_ϕ[Quantum]
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus
    loss =
        (1 // 2) *
        im *
        (
            bar(c) * bar(q) * (c(minus)^2 + q(minus)^2) -
            c(plus) * q(plus) * (bar(c)^2 + bar(q)^2) +
            2 * bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
        )
    L = InteractionLagrangian(loss, :γ)
    G = DressedPropagator(L, Val(1), Val(3); simplify=true, _set_reg_to_zero=false)
    GF = fourier_transform(G)
    ΣW = wigner_transform(SelfEnergy(GF); gradient_order=Val(0))
    spectral = spectral_dispersive_collision(
        off_shell_collision_expression(kinetic_expression(ΣW))
    )
    return KC.canonical_frequency_collision(spectral)
end

@testset "generated first-order loss collision reduction" begin
    canonical = generated_collision_reduction_fixture()
    reduced = @inferred KC.reduce_frequency_collision(canonical)
    @test isempty(KC.reduced_blocked_terms(reduced))
    @test isempty(KC.reduced_trotter_terms(reduced))
    @test length(KC.reduced_regular_terms(reduced)) == 1

    statistical = only(values(KC.reduced_regular_terms(reduced)))
    occupation = @inferred KC.occupation_reduced_expression(reduced)
    @test length(KC.occupation_reduced_terms(occupation)) == 1
    occupation_polynomial = only(values(KC.occupation_reduced_terms(occupation)))
    @test length(occupation_polynomial) == 1
    monomial, coefficient = only(occupation_polynomial)
    @test length(monomial) == 2
    @test coefficient == -4 // 1

    @test length(statistical) == 4
end
