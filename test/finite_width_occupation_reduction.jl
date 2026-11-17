using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields finite_width_occupation_ϕ::Boson

function generated_finite_width_occupation_problem()
    c = finite_width_occupation_ϕ[Classical]
    q = finite_width_occupation_ϕ[Quantum]
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

function finite_width_occupation_lines(collision)
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

function finite_width_unit_model(collision)
    return KC.LorentzianSpectralModel(
        Dict(
            line => KC.LorentzianSpectralData(0 // 1, 1 // 1) for
            line in finite_width_occupation_lines(collision)
        ),
    )
end

function expected_finite_width_occupation(evaluated)
    S = KC.statistics(evaluated)
    target = KC.target_family(evaluated)
    out = Dict()

    for (expression, include_external) in (
        KC.collision_offset(evaluated) => false,
        KC.collision_distribution_coefficient(evaluated) => true,
    )
        for (term, coefficient) in KC.finite_width_resolved_terms(expression)
            statistical = KC._collision_statistical_monomial(
                KC.finite_width_source_term(term), target, include_external
            )
            C = typeof(coefficient)
            polynomial = KC.StatisticalPolynomial{C,S}([statistical => coefficient])
            occupation = KC.occupation_collision_polynomial(polynomial)
            if haskey(out, term)
                out[term] = out[term] + occupation
                iszero(out[term]) && delete!(out, term)
            elseif !iszero(occupation)
                out[term] = occupation
            end
        end
    end
    return out
end

@testset "resolved finite-width collision lowers with common occupation algebra" begin
    frequency_collision = generated_finite_width_occupation_problem()
    model = finite_width_unit_model(frequency_collision)
    evaluated = @inferred KC.evaluate_finite_width_collision(
        frequency_collision, model, 0 // 1
    )
    occupation = @inferred KC.finite_width_occupation_collision(evaluated)

    @test occupation isa KC.FiniteWidthOccupationCollision
    @test KC.order(occupation) == KC.order(evaluated)
    @test KC.statistics(occupation) === Boson
    @test KC.parameters(occupation) == KC.parameters(evaluated)
    @test KC.target_family(occupation) === finite_width_occupation_ϕ
    @test KC.gradient_order(occupation) == Val(0)
    @test KC.wigner_context(occupation) == KC.wigner_context(evaluated)

    @test !isempty(KC.finite_width_occupation_terms(occupation))
    @test KC.finite_width_occupation_terms(occupation) ==
        expected_finite_width_occupation(evaluated)

    @test KC.finite_width_unsupported_offset_terms(occupation) ==
        KC.finite_width_unsupported_terms(KC.collision_offset(evaluated))
    @test KC.finite_width_unsupported_distribution_terms(occupation) ==
        KC.finite_width_unsupported_terms(
        KC.collision_distribution_coefficient(evaluated)
    )
end
