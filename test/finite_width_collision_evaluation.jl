using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields finite_width_eval_ϕ::Boson

function generated_finite_width_evaluation_problem()
    c = finite_width_eval_ϕ[Classical]
    q = finite_width_eval_ϕ[Quantum]
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
    spectral = @inferred spectral_dispersive_collision(offshell)
    return @inferred KC.finite_width_frequency_collision(spectral)
end

function resolved_lines(collision)
    lines = KC.SpectralLineIdentity{Boson}[]
    for expression in
        (KC.collision_offset(collision), KC.collision_distribution_coefficient(collision))
        for (term, coefficient) in expression
            iszero(coefficient) && continue
            kind = KC.finite_width_reduction_kind(term)
            if kind === KC.FiniteWidthFactorized
                for factor in KC.spectral_factors(KC.finite_width_spectral_reduction(term))
                    push!(lines, KC.spectral_line(factor))
                end
            elseif kind === KC.FiniteWidthConvolution
                reduction = KC.finite_width_convolution_reduction(term)
                append!(lines, KC.convolution_pivot_lines(reduction))
                push!(lines, KC.convolution_dependent_line(reduction))
            end
        end
    end
    unique!(lines)
    return lines
end

function unit_width_model(collision)
    data = Dict(
        line => KC.LorentzianSpectralData(0 // 1, 1 // 1) for
        line in resolved_lines(collision)
    )
    return KC.LorentzianSpectralModel(data)
end

function expected_resolved(expression, model, ω_external)
    out = Dict()
    for (term, coefficient) in expression
        kind = KC.finite_width_reduction_kind(term)
        if kind === KC.FiniteWidthFactorized
            weight = KC.evaluate_spectral_weight(
                KC.finite_width_spectral_reduction(term), model
            )
            out[term] = coefficient * weight
        elseif kind === KC.FiniteWidthConvolution
            weight = KC.evaluate_spectral_convolution(
                KC.finite_width_convolution_reduction(term), model, ω_external
            )
            out[term] = coefficient * weight
        end
    end
    return out
end

@testset "finite-width collision evaluates resolved frequency weights" begin
    collision = generated_finite_width_evaluation_problem()
    model = unit_width_model(collision)
    ω_external = 0 // 1
    evaluated = @inferred KC.evaluate_finite_width_collision(collision, model, ω_external)

    @test evaluated isa KC.FiniteWidthEvaluatedCollision
    @test KC.order(evaluated) == KC.order(collision)
    @test KC.statistics(evaluated) === Boson
    @test KC.parameters(evaluated) == KC.parameters(collision)
    @test KC.target_family(evaluated) === finite_width_eval_ϕ
    @test KC.gradient_order(evaluated) == Val(0)
    @test KC.wigner_context(evaluated) == KC.wigner_context(collision)

    for (source, result) in (
        KC.collision_offset(collision) => KC.collision_offset(evaluated),
        KC.collision_distribution_coefficient(collision) =>
            KC.collision_distribution_coefficient(evaluated),
    )
        @test KC.finite_width_resolved_terms(result) ==
            expected_resolved(source, model, ω_external)
        for (term, coefficient) in source
            if KC.finite_width_reduction_kind(term) === KC.FiniteWidthUnsupported
                @test KC.finite_width_unsupported_terms(result)[term] == coefficient
            end
        end
    end

    resolved = [
        term for expression in
        (KC.collision_offset(collision), KC.collision_distribution_coefficient(collision))
        for (term, coefficient) in expression if !iszero(coefficient) &&
            KC.finite_width_reduction_kind(term) === KC.FiniteWidthFactorized &&
            any(
                factor -> KC.spectral_multiplicity(factor) == 2,
                KC.spectral_factors(KC.finite_width_spectral_reduction(term)),
            )
    ]
    @test !isempty(resolved)
    repeated_term = first(resolved)
    repeated_line = KC.spectral_line(
        only(
            factor for factor in
            KC.spectral_factors(KC.finite_width_spectral_reduction(repeated_term)) if
            KC.spectral_multiplicity(factor) == 2
        ),
    )

    doubled_data = Dict(
        line => KC.LorentzianSpectralData(0 // 1, line == repeated_line ? 2 // 1 : 1 // 1)
        for line in resolved_lines(collision)
    )
    doubled_model = KC.LorentzianSpectralModel(doubled_data)
    doubled = @inferred KC.evaluate_finite_width_collision(
        collision, doubled_model, ω_external
    )

    function resolved_value(result, term)
        offset = KC.finite_width_resolved_terms(KC.collision_offset(result))
        haskey(offset, term) && return offset[term]
        return KC.finite_width_resolved_terms(
            KC.collision_distribution_coefficient(result)
        )[term]
    end

    @test 2 * resolved_value(doubled, repeated_term) ==
        resolved_value(evaluated, repeated_term)

    empty_model = KC.LorentzianSpectralModel(
        Dict{KC.SpectralLineIdentity{Boson},KC.LorentzianSpectralData{Rational{Int64}}}()
    )
    @test_throws ArgumentError KC.evaluate_finite_width_collision(
        collision, empty_model, ω_external
    )
end

@testset "finite-width evaluated collision public boundary" begin
    for name in (
        :FiniteWidthEvaluatedExpression,
        :FiniteWidthEvaluatedCollision,
        :finite_width_resolved_terms,
        :finite_width_unsupported_terms,
        :evaluate_finite_width_collision,
    )
        if VERSION >= v"1.11"
            @test !Base.isexported(KC, name)
            @test Base.ispublic(KC, name)
        else
            @test !(name in Base.names(KC; all=false, imported=false))
        end
    end
end
