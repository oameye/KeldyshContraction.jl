using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields finite_width_collision_ϕ::Boson

function generated_second_order_finite_width_collision()
    c = finite_width_collision_ϕ[Classical]
    q = finite_width_collision_ϕ[Quantum]
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
    offshell = @inferred off_shell_collision_expression(KΣ)
    return @inferred spectral_dispersive_collision(offshell)
end

function source_terms(expression)
    return Dict(
        KC.finite_width_source_term(term) => coefficient for
        (term, coefficient) in expression
    )
end

function all_finite_width_terms(collision)
    return [
        term for expression in
        (KC.collision_offset(collision), KC.collision_distribution_coefficient(collision))
        for (term, coefficient) in expression if !iszero(coefficient)
    ]
end

@testset "whole-collision finite-width frequency classification" begin
    spectral = generated_second_order_finite_width_collision()
    collision = @inferred KC.finite_width_frequency_collision(spectral)

    @test collision isa KC.FiniteWidthFrequencyCollision
    @test KC.order(collision) == KC.order(spectral)
    @test KC.statistics(collision) === Boson
    @test KC.parameters(collision) == KC.parameters(spectral)
    @test KC.target_family(collision) === finite_width_collision_ϕ
    @test KC.gradient_order(collision) == Val(0)
    @test KC.wigner_context(collision) == KC.wigner_context(spectral)

    # Classification is lossless: every exact source term and coefficient remains in the same
    # affine Kadanoff--Baym component.
    @test source_terms(KC.collision_offset(collision)) ==
        Dict(KC.collision_offset(spectral))
    @test source_terms(KC.collision_distribution_coefficient(collision)) ==
        Dict(KC.collision_distribution_coefficient(spectral))

    terms = all_finite_width_terms(collision)
    factorized = [
        term for
        term in terms if KC.finite_width_reduction_kind(term) === KC.FiniteWidthFactorized
    ]
    convolutions = [
        term for
        term in terms if KC.finite_width_reduction_kind(term) === KC.FiniteWidthConvolution
    ]
    unsupported = [
        term for
        term in terms if KC.finite_width_reduction_kind(term) === KC.FiniteWidthUnsupported
    ]

    @test !isempty(factorized)
    @test !isempty(convolutions)
    @test !isempty(unsupported)
    @test all(term -> KC.target_family(term) === finite_width_collision_ϕ, terms)

    # Eq. 55b: a generated repeated spectral line is resolved by physical line identity and
    # retains the explicit γ² / Γ enhancement.
    repeated = [
        term for term in factorized if any(
            factor -> KC.spectral_multiplicity(factor) >= 2,
            KC.spectral_factors(KC.finite_width_spectral_reduction(term)),
        )
    ]
    @test !isempty(repeated)
    repeated_term = first(repeated)
    repeated_reduction = KC.finite_width_spectral_reduction(repeated_term)
    @test KC.spectral_reduction_resolved(repeated_reduction)
    @test any(
        factor ->
            KC.spectral_multiplicity(factor) == 2 && KC.linewidth_exponent(factor) == -1,
        KC.spectral_factors(repeated_reduction),
    )
    counting = KC.width_aware_power_counting(KC.parameters(collision), repeated_reduction)
    @test KC.parameters(counting) == KC.ParameterMonomial(:γ)^2
    @test any(power -> KC.linewidth_exponent(power) == -1, KC.linewidth_powers(counting))

    # Eq. 55a: the generated regular sector is the existing one-residual-shell Cauchy
    # convolution, with exactly the same mismatch and Jacobian as the strict-QP shell limit.
    convolution_term = first(convolutions)
    convolution = KC.finite_width_convolution_reduction(convolution_term)
    @test KC.convolution_reduction_resolved(convolution)
    source = KC.finite_width_source_term(convolution_term)
    strict = KC.full_rank_frequency_reduction(source, finite_width_collision_ϕ)
    shell, support_factor = KC.energy_shell(KC.convolution_energy_mismatch(convolution))
    support = KC.frequency_support(strict)
    @test support.shells == [shell]
    @test isempty(support.principal_values)
    @test KC.frequency_factor(strict) ==
        KC.convolution_jacobian(convolution) * support_factor

    # Unsupported structures remain concrete source terms rather than acquiring an implicit
    # cutoff, sharp-shell value, or fallback finite number.
    unsupported_term = first(unsupported)
    @test KC.finite_width_source_term(unsupported_term) isa KC.SpectralDispersiveTerm
    @test_throws ArgumentError KC.finite_width_spectral_reduction(unsupported_term)
    @test_throws ArgumentError KC.finite_width_convolution_reduction(unsupported_term)
end

@testset "finite-width frequency collision public boundary" begin
    names = (
        :FiniteWidthFrequencyReductionKind,
        :FiniteWidthFactorized,
        :FiniteWidthConvolution,
        :FiniteWidthUnsupported,
        :FiniteWidthFrequencyTerm,
        :FiniteWidthFrequencyExpression,
        :FiniteWidthFrequencyCollision,
        :finite_width_source_term,
        :finite_width_reduction_kind,
        :finite_width_spectral_reduction,
        :finite_width_convolution_reduction,
        :finite_width_frequency_collision,
    )
    for name in names
        if VERSION >= v"1.11"
            @test !Base.isexported(KC, name)
            @test Base.ispublic(KC, name)
        else
            @test !(name in Base.names(KC; all=false, imported=false))
        end
    end
end
