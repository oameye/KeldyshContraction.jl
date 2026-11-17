using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields finite_width_response_ϕ::Boson

function generated_finite_width_response_problem()
    c = finite_width_response_ϕ[Classical]
    q = finite_width_response_ϕ[Quantum]
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
    G = DressedPropagator(L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false)
    ΣF = SelfEnergy(fourier_transform(G))
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    kinetic = kinetic_expression(ΣW)
    offshell = off_shell_collision_expression(kinetic)
    spectral = spectral_dispersive_collision(offshell)
    return KC.finite_width_frequency_collision(spectral)
end

function finite_width_response_lines(collision)
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

function finite_width_response_model(collision)
    return KC.LorentzianSpectralModel(
        Dict(
            line => KC.LorentzianSpectralData(0 // 1, 1 // 1, 1 // 1) for
            line in finite_width_response_lines(collision)
        ),
    )
end

function total_occupation_polynomial(collision, selected)
    pieces = Any[]
    target = KC.target_family(collision)
    for (expression, include_external) in (
        KC.collision_offset(collision) => false,
        KC.collision_distribution_coefficient(collision) => true,
    )
        for (term, coefficient) in expression
            isequal(term, selected) || continue
            push!(
                pieces,
                KC._finite_width_occupation_polynomial(
                    term, coefficient, target, include_external
                ),
            )
        end
    end
    isempty(pieces) && error("selected finite-width term is absent from collision")
    return reduce(+, pieces)
end

function fixed_model_background_response(collision, model, background)
    evaluated = KC.evaluate_finite_width_collision(collision, model, 0 // 1)
    occupation = KC.finite_width_occupation_collision(evaluated)
    out = Dict()
    for (term, polynomial) in KC.finite_width_occupation_terms(occupation)
        linearization = KC.evaluate_occupation_linearization(
            KC.occupation_linearization(polynomial), background
        )
        iszero(linearization) || (out[term] = linearization)
    end
    return out, occupation
end

function repeated_response_term(collision, background)
    seen = Set{Any}()
    for expression in
        (KC.collision_offset(collision), KC.collision_distribution_coefficient(collision))
        for (term, _) in expression
            term in seen && continue
            push!(seen, term)
            KC.finite_width_reduction_kind(term) === KC.FiniteWidthFactorized || continue
            reduction = KC.finite_width_spectral_reduction(term)
            factors = KC.spectral_factors(reduction)
            count(factor -> KC.spectral_multiplicity(factor) == 2, factors) == 1 || continue
            all(factor -> KC.spectral_multiplicity(factor) in (1, 2), factors) || continue
            polynomial = total_occupation_polynomial(collision, term)
            iszero(KC.evaluate_occupation_polynomial(polynomial, background)) && continue
            linearization = KC.occupation_linearization(polynomial)
            iszero(linearization) && continue
            return term, polynomial, linearization
        end
    end
    return error("generated collision contains no usable repeated-line response term")
end

@testset "zero spectral response reproduces fixed-model linearization" begin
    collision = generated_finite_width_response_problem()
    model = finite_width_response_model(collision)
    background = KC.OccupationBackground(_ -> 1 // 1, Rational{Int})
    Variation = KC.LorentzianSpectralModelVariation{Rational{Int},Boson}
    empty_response = KC.OccupationSpectralResponse(
        Dict{KC.OccupationAtom{Boson},Variation}()
    )

    response = @inferred KC.linearize_finite_width_collision(
        collision, model, background, empty_response, 0 // 1
    )
    expected, occupation = fixed_model_background_response(collision, model, background)

    @test response isa KC.FiniteWidthCollisionLinearization
    @test KC.order(response) == KC.order(collision)
    @test KC.statistics(response) === Boson
    @test KC.parameters(response) == KC.parameters(collision)
    @test KC.target_family(response) === finite_width_response_ϕ
    @test KC.gradient_order(response) == Val(0)
    @test KC.wigner_context(response) == KC.wigner_context(collision)
    @test isempty(KC.finite_width_spectral_response_terms(response))

    combined = KC.finite_width_linearized_terms(response)
    occupation_response = KC.finite_width_occupation_response_terms(response)
    @test Set(keys(combined)) == Set(keys(expected))
    @test Set(keys(occupation_response)) == Set(keys(expected))
    for term in keys(expected)
        expected_terms = KC.background_occupation_linearization_terms(expected[term])
        @test KC.background_occupation_linearization_terms(combined[term]) == expected_terms
        @test KC.background_occupation_linearization_terms(occupation_response[term]) ==
            expected_terms
    end

    @test KC.finite_width_unsupported_offset_terms(response) ==
        KC.finite_width_unsupported_offset_terms(occupation)
    @test KC.finite_width_unsupported_distribution_terms(response) ==
        KC.finite_width_unsupported_distribution_terms(occupation)
end

@testset "generated repeated line obeys exact product rule" begin
    collision = generated_finite_width_response_problem()
    model = finite_width_response_model(collision)
    background = KC.OccupationBackground(_ -> 1 // 1, Rational{Int})
    term, polynomial, linearization = repeated_response_term(collision, background)
    reduction = KC.finite_width_spectral_reduction(term)
    repeated = only(
        factor for
        factor in KC.spectral_factors(reduction) if KC.spectral_multiplicity(factor) == 2
    )
    repeated_line = KC.spectral_line(repeated)

    background_linearization = KC.evaluate_occupation_linearization(
        linearization, background
    )
    atom = first(
        first(KC.background_occupation_linearization_terms(background_linearization))
    )

    lines = finite_width_response_lines(collision)
    variations = Dict(
        line => KC.LorentzianSpectralDataVariation(
            0 // 1, isequal(line, repeated_line) ? 1 // 1 : 0 // 1, 0 // 1
        ) for line in lines
    )
    model_variation = KC.LorentzianSpectralModelVariation(variations)
    spectral_response = KC.OccupationSpectralResponse(Dict(atom => model_variation))
    response = @inferred KC.linearize_finite_width_collision(
        collision, model, background, spectral_response, 0 // 1
    )

    jacobian = KC.spectral_jacobian(reduction)
    weight = 2 * jacobian
    δweight = -2 * jacobian
    background_polynomial = KC.evaluate_occupation_polynomial(polynomial, background)
    derivative_polynomial = Dict(
        KC.background_occupation_linearization_terms(background_linearization)
    )[atom]

    occupation_terms = Dict(
        KC.background_occupation_linearization_terms(
            KC.finite_width_occupation_response_terms(response)[term]
        ),
    )
    spectral_terms = Dict(
        KC.background_occupation_linearization_terms(
            KC.finite_width_spectral_response_terms(response)[term]
        ),
    )
    combined_terms = Dict(
        KC.background_occupation_linearization_terms(
            KC.finite_width_linearized_terms(response)[term]
        ),
    )

    @test occupation_terms[atom] == weight * derivative_polynomial
    @test spectral_terms[atom] == background_polynomial * δweight
    @test combined_terms[atom] == occupation_terms[atom] + spectral_terms[atom]
end

@testset "finite-width response API is qualified-public" begin
    for name in (
        :OccupationBackground,
        :OccupationLinearization,
        :BackgroundOccupationLinearization,
        :OccupationSpectralResponse,
        :FiniteWidthCollisionLinearization,
        :occupation_background_value,
        :evaluate_occupation_polynomial,
        :occupation_linearization,
        :evaluate_occupation_linearization,
        :spectral_occupation_response_terms,
        :finite_width_linearized_terms,
        :finite_width_occupation_response_terms,
        :finite_width_spectral_response_terms,
        :linearize_finite_width_collision,
    )
        if VERSION >= v"1.11"
            @test !Base.isexported(KC, name)
            @test Base.ispublic(KC, name)
        else
            @test !(name in Base.names(KC; all=false, imported=false))
        end
    end
end
