using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields generated_width_response_ϕ::Boson

function generated_width_response_loss_lagrangian()
    c = generated_width_response_ϕ[Classical]
    q = generated_width_response_ϕ[Quantum]
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

function generated_width_response_kinetic(
    order::Val{O}, topological_order::Val{T}
) where {O,T}
    L = generated_width_response_loss_lagrangian()
    G = DressedPropagator(
        L, order, topological_order; simplify=true, _set_reg_to_zero=false
    )
    ΣF = SelfEnergy(fourier_transform(G))
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    return kinetic_expression(ΣW)
end

function generated_width_response_linewidth_kernel()
    kinetic = generated_width_response_kinetic(Val(1), Val(3))
    return @inferred KC.spectral_self_energy_kernel(kinetic)
end

function generated_width_response_collision()
    kinetic = generated_width_response_kinetic(Val(2), Val(5))
    offshell = off_shell_collision_expression(kinetic)
    spectral = spectral_dispersive_collision(offshell)
    return KC.finite_width_frequency_collision(spectral)
end

function generated_width_response_lines(collision)
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

function generated_width_response_polynomial(collision, selected)
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

function generated_width_repeated_term(collision, background)
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
            polynomial = generated_width_response_polynomial(collision, term)
            iszero(KC.evaluate_occupation_polynomial(polynomial, background)) && continue
            linearization = KC.evaluate_occupation_linearization(
                KC.occupation_linearization(polynomial), background
            )
            iszero(linearization) && continue
            return term, polynomial, linearization
        end
    end
    return error("generated collision contains no usable repeated-line response term")
end

function generated_scalar_linewidth(kernel, occupation_value, γ_value)
    background = KC.OccupationBackground(_ -> occupation_value, Rational{Int})
    evaluated = KC.evaluate_spectral_self_energy_background(kernel, background)
    sector, coefficient = only(KC.background_spectral_self_energy_terms(evaluated))

    # This benchmark closure is deliberately local: KC still exposes the routed loop integral.
    # Here ∫_q 1 = 1, so the generated 4 n_q source gives Γ = 4γ n and δΓ = 4γ δn.
    unit_normalized_loop_integral = (_, value) -> real(value)
    return γ_value * unit_normalized_loop_integral(sector, coefficient)
end

@testset "generated linewidth closes self-consistent finite-width loss response" begin
    linewidth_kernel = generated_width_response_linewidth_kernel()
    collision = generated_width_response_collision()

    @test KC.parameters(linewidth_kernel) == KC.ParameterMonomial(:γ)
    @test isempty(KC.spectral_self_energy_blocked_terms(linewidth_kernel))
    @test isempty(KC.spectral_self_energy_trotter_terms(linewidth_kernel))

    γ_value = 1 // 4
    background_value = 3 // 1
    unit_variation = 1 // 1
    Γbar = generated_scalar_linewidth(linewidth_kernel, background_value, γ_value)
    δΓ = generated_scalar_linewidth(linewidth_kernel, unit_variation, γ_value)
    @test Γbar == 3 // 1
    @test δΓ == 1 // 1

    background = KC.OccupationBackground(_ -> background_value, Rational{Int})
    term, polynomial, background_linearization = generated_width_repeated_term(
        collision, background
    )
    atom = first(
        first(KC.background_occupation_linearization_terms(background_linearization))
    )

    lines = generated_width_response_lines(collision)
    model = KC.LorentzianSpectralModel(
        Dict(line => KC.LorentzianSpectralData(0 // 1, Γbar, 1 // 1) for line in lines)
    )
    model_variation = KC.LorentzianSpectralModelVariation(
        Dict(
            line => KC.LorentzianSpectralDataVariation(0 // 1, δΓ, 0 // 1) for line in lines
        ),
    )
    Variation = KC.LorentzianSpectralModelVariation{Rational{Int},Boson}
    frozen_spectral_response = KC.OccupationSpectralResponse(
        Dict{KC.OccupationAtom{Boson},Variation}()
    )
    full_spectral_response = KC.OccupationSpectralResponse(Dict(atom => model_variation))

    frozen = @inferred KC.linearize_finite_width_collision(
        collision, model, background, frozen_spectral_response, 0 // 1
    )
    full = @inferred KC.linearize_finite_width_collision(
        collision, model, background, full_spectral_response, 0 // 1
    )

    reduction = KC.finite_width_spectral_reduction(term)
    repeated = only(
        factor for
        factor in KC.spectral_factors(reduction) if KC.spectral_multiplicity(factor) == 2
    )
    @test KC.spectral_multiplicity(repeated) == 2

    jacobian = KC.spectral_jacobian(reduction)
    expected_weight = 2 * jacobian / Γbar
    expected_δweight = -2 * jacobian * δΓ / Γbar^2
    Pbar = KC.evaluate_occupation_polynomial(polynomial, background)
    δP = Dict(KC.background_occupation_linearization_terms(background_linearization))[atom]

    frozen_terms = Dict(
        KC.background_occupation_linearization_terms(
            KC.finite_width_linearized_terms(frozen)[term]
        ),
    )
    occupation_terms = Dict(
        KC.background_occupation_linearization_terms(
            KC.finite_width_occupation_response_terms(full)[term]
        ),
    )
    spectral_terms = Dict(
        KC.background_occupation_linearization_terms(
            KC.finite_width_spectral_response_terms(full)[term]
        ),
    )
    full_terms = Dict(
        KC.background_occupation_linearization_terms(
            KC.finite_width_linearized_terms(full)[term]
        ),
    )

    @test frozen_terms[atom] == expected_weight * δP
    @test occupation_terms[atom] == frozen_terms[atom]
    @test spectral_terms[atom] == Pbar * expected_δweight
    @test full_terms[atom] == frozen_terms[atom] + spectral_terms[atom]
    @test full_terms[atom] != frozen_terms[atom]

    @test KC.finite_width_unsupported_offset_terms(full) ==
        KC.finite_width_unsupported_offset_terms(frozen)
    @test KC.finite_width_unsupported_distribution_terms(full) ==
        KC.finite_width_unsupported_distribution_terms(frozen)
end
