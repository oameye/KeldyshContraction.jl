using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields offshell_loss2_ϕ::Boson

function generated_second_order_off_shell_collision()
    c = offshell_loss2_ϕ[Classical]
    q = offshell_loss2_ϕ[Quantum]
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

    L = @inferred InteractionLagrangian(loss, :γ)
    G = @inferred DressedPropagator(
        L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false
    )
    ΣF = @inferred SelfEnergy(fourier_transform(G))
    ΣW = @inferred wigner_transform(ΣF; gradient_order=Val(0))
    KΣ = @inferred kinetic_expression(ΣW)
    return @inferred off_shell_collision_expression(KΣ)
end

function spectral_momenta(term)
    return [
        KC.momentum(line) for
        (line, kind) in zip(KC.kinetic_lines(term.carrier), KC.spectral_dispersive_kinds(term)) if
        kind === KC.CollisionSpectral
    ]
end

function maximum_spectral_multiplicity(term)
    momenta = spectral_momenta(term)
    isempty(momenta) && return 0
    return maximum(count(==(momentum), momenta) for momentum in momenta)
end

@testset "generated γ² off-shell loss keeps regular and repeated spectral sectors" begin
    collision = generated_second_order_off_shell_collision()

    @test collision isa KC.OffShellCollisionExpression
    @test @inferred(KC.order(collision)) == 2
    @test @inferred(KC.statistics(collision)) === Boson
    @test @inferred(KC.parameters(collision)) == KC.ParameterMonomial(:γ)^2
    @test @inferred(KC.target_family(collision)) === offshell_loss2_ϕ
    @test @inferred(KC.gradient_order(collision)) == Val(0)
    @test !iszero(@inferred(KC.collision_offset(collision)))
    @test !iszero(@inferred(KC.collision_distribution_coefficient(collision)))

    # This is an exact basis change of the same off-shell collision functional. No shell delta
    # functions, occupation substitution, or finite-width prescription has been introduced.
    spectral = @inferred spectral_dispersive_collision(collision)
    terms = [
        term for expression in (
            KC.collision_offset(spectral), KC.collision_distribution_coefficient(spectral)
        ) for (term, coefficient) in expression if !iszero(coefficient) &&
        count(==(KC.CollisionSpectral), KC.spectral_dispersive_kinds(term)) == 3
    ]

    # The ordinary sunset contribution has three distinct spectral momenta (Eq. 55a), while
    # the repeated-line contribution contains A(q)^2 A(q') already off shell (Eq. 55b).
    @test any(term -> maximum_spectral_multiplicity(term) == 1, terms)
    @test any(term -> maximum_spectral_multiplicity(term) >= 2, terms)

    # Only the strict quasiparticle reduction turns the repeated spectral line into a blocker.
    # The regular sector remains finite and the existing reduction semantics are unchanged.
    reduced = @inferred reduce_frequency_collision(spectral)
    @test length(KC.reduced_regular_terms(reduced)) == 1
    @test !isempty(KC.reduced_blocked_terms(reduced))
    @test isempty(KC.reduced_trotter_terms(reduced))
end
