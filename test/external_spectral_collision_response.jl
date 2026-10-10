using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields external_response_ϕ::Boson

function external_response_loss_lagrangian()
    c = external_response_ϕ[Classical]
    q = external_response_ϕ[Quantum]
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
    return InteractionLagrangian(loss, :γ)
end

function generated_external_response_collision()
    L = external_response_loss_lagrangian()
    G = DressedPropagator(L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false)
    ΣF = SelfEnergy(fourier_transform(G))
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    kinetic = kinetic_expression(ΣW)
    offshell = off_shell_collision_expression(kinetic)
    spectral = spectral_dispersive_collision(offshell)
    return KC.finite_width_frequency_collision(spectral)
end

function generated_external_response_convolution(collision)
    for (include_external, expression) in (
        false => KC.collision_offset(collision),
        true => KC.collision_distribution_coefficient(collision),
    )
        for (term, coefficient) in expression
            iszero(coefficient) && continue
            KC.finite_width_reduction_kind(term) === KC.FiniteWidthConvolution || continue
            return term, coefficient, include_external
        end
    end
    return error("generated collision contains no resolved convolution sector")
end

function one_term_external_response_collision(
    collision, term, coefficient, include_external
)
    offset = KC.collision_offset(collision)
    distribution = KC.collision_distribution_coefficient(collision)
    Term = typeof(term)
    Coefficient = typeof(coefficient)
    offset_terms = Dict{Term,Coefficient}()
    distribution_terms = Dict{Term,Coefficient}()
    if include_external
        distribution_terms[term] = coefficient
    else
        offset_terms[term] = coefficient
    end
    offset_one = typeof(offset)(offset_terms, KC.wigner_context(offset))
    distribution_one = typeof(distribution)(
        distribution_terms, KC.wigner_context(distribution)
    )
    return typeof(collision)(
        offset_one,
        distribution_one,
        KC.target_family(collision),
        KC.parameters(collision),
        KC.wigner_context(collision),
    )
end

@testset "generated sunset distinguishes fixed-ω and external spectral response" begin
    collision = generated_external_response_collision()
    term, coefficient, include_external = generated_external_response_convolution(collision)
    collision = one_term_external_response_collision(
        collision, term, coefficient, include_external
    )
    target = KC.target_family(collision)
    reduction = KC.finite_width_convolution_reduction(term)
    external_line = KC._finite_width_external_spectral_line(term, target)
    source = KC.finite_width_source_term(term)
    source_basis = KC.momentum_basis(source)
    external_index = KC._external_frequency_basis_index(source)
    internal_lines = vcat(
        KC.convolution_pivot_lines(reduction), [KC.convolution_dependent_line(reduction)]
    )
    lines = unique(vcat(internal_lines, [external_line]))

    @test source_basis[external_index] == KC.external_wigner_momentum(source)
    @test KC.spectral_line_momentum(external_line) ==
        KC.basis_momentum(source_basis, external_index)
    @test length(internal_lines) == 3
    @test length(lines) == 4
    @test abs(KC.convolution_external_coefficient(reduction)) == 1 // 1
    @test all(
        coefficient -> abs(coefficient) == 1 // 1, KC.convolution_coefficients(reduction)
    )

    model = KC.LorentzianSpectralModel(
        Dict(line => KC.LorentzianSpectralData(0 // 1, 1 // 1, 1 // 1) for line in lines)
    )
    variation = KC.LorentzianSpectralModelVariation(
        Dict(
            line => KC.LorentzianSpectralDataVariation(0 // 1, 1 // 1, 0 // 1) for
            line in lines
        ),
    )
    background = KC.OccupationBackground(_ -> 3 // 1, Rational{Int})
    polynomial = KC._finite_width_occupation_polynomial(
        term, coefficient, target, include_external
    )
    background_linearization = KC.evaluate_occupation_linearization(
        KC.occupation_linearization(polynomial), background
    )
    atom, δP = first(KC.background_occupation_linearization_terms(background_linearization))
    Pbar = KC.evaluate_occupation_polynomial(polynomial, background)
    response = KC.OccupationSpectralResponse(Dict(atom => variation))

    externally_projected = @inferred KC.linearize_external_spectral_collision(
        collision, model, background, response
    )
    fixed_external_frequency = @inferred KC.linearize_finite_width_collision(
        collision, model, background, response, 0 // 1
    )

    J = KC.convolution_jacobian(reduction)
    @test KC.convolution_effective_linewidth(reduction, model) == 3 // 1
    @test KC.external_spectral_projection_linewidth(reduction, model, external_line) ==
        4 // 1
    @test KC.external_spectral_projection_mismatch(reduction, model, external_line) ==
        0 // 1
    @test KC.evaluate_spectral_convolution(reduction, model, 0 // 1) == 4J / 3
    @test KC.evaluate_external_spectral_projection(reduction, model, external_line) == J
    @test KC.evaluate_spectral_convolution_variation(reduction, model, variation, 0 // 1) ==
        -4J / 3
    @test KC.evaluate_external_spectral_projection_variation(
        reduction, model, variation, external_line
    ) == -J

    external_occupation = Dict(
        KC.background_occupation_linearization_terms(
            KC.external_spectral_occupation_response_terms(externally_projected)[term]
        ),
    )
    external_spectral = Dict(
        KC.background_occupation_linearization_terms(
            KC.external_spectral_model_response_terms(externally_projected)[term]
        ),
    )
    external_full = Dict(
        KC.background_occupation_linearization_terms(
            KC.external_spectral_linearized_terms(externally_projected)[term]
        ),
    )
    fixed_occupation = Dict(
        KC.background_occupation_linearization_terms(
            KC.finite_width_occupation_response_terms(fixed_external_frequency)[term]
        ),
    )
    fixed_spectral = Dict(
        KC.background_occupation_linearization_terms(
            KC.finite_width_spectral_response_terms(fixed_external_frequency)[term]
        ),
    )

    @test external_occupation[atom] == J * δP
    @test external_spectral[atom] == -J * Pbar
    @test external_full[atom] == J * δP - J * Pbar
    @test fixed_occupation[atom] == (4J / 3) * δP
    @test fixed_spectral[atom] == (-4J / 3) * Pbar
    @test external_full[atom] != Dict(
        KC.background_occupation_linearization_terms(
            KC.finite_width_linearized_terms(fixed_external_frequency)[term]
        ),
    )[atom]

    missing_external = KC.LorentzianSpectralModelVariation(
        Dict(
            line => KC.LorentzianSpectralDataVariation(0 // 1, 1 // 1, 0 // 1) for
            line in internal_lines
        ),
    )
    missing_response = KC.OccupationSpectralResponse(Dict(atom => missing_external))
    @test_throws ArgumentError KC.linearize_external_spectral_collision(
        collision, model, background, missing_response
    )
end
