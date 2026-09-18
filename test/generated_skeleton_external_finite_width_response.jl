using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields skeleton_external_response_ϕ::Boson

function skeleton_external_response_loss_lagrangian()
    c = skeleton_external_response_ϕ[Classical]
    q = skeleton_external_response_ϕ[Quantum]
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

function skeleton_external_response_kinetic(
    order::Val{O}, topological_order::Val{T}; skeleton
) where {O,T}
    L = skeleton_external_response_loss_lagrangian()
    G = DressedPropagator(
        L, order, topological_order; simplify=true, _set_reg_to_zero=false
    )
    GF = fourier_transform(G)
    ΣF = skeleton ? KC.skeleton_self_energy(GF) : SelfEnergy(GF)
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    return kinetic_expression(ΣW)
end

function skeleton_external_response_collision(; skeleton)
    kinetic = skeleton_external_response_kinetic(Val(2), Val(5); skeleton)
    offshell = off_shell_collision_expression(kinetic)
    spectral = spectral_dispersive_collision(offshell)
    return KC.finite_width_frequency_collision(spectral)
end

function skeleton_external_response_linewidth_kernel()
    kinetic = skeleton_external_response_kinetic(Val(1), Val(3); skeleton=true)
    return KC.spectral_self_energy_kernel(kinetic)
end

function skeleton_external_response_kinds(collision)
    kinds = typeof(KC.FiniteWidthFactorized)[]
    for expression in
        (KC.collision_offset(collision), KC.collision_distribution_coefficient(collision))
        for (term, coefficient) in expression
            iszero(coefficient) && continue
            push!(kinds, KC.finite_width_reduction_kind(term))
        end
    end
    return kinds
end

function skeleton_external_response_polynomial(collision, selected)
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

function skeleton_external_response_convolution(collision, background)
    seen = Set{Any}()
    for expression in
        (KC.collision_offset(collision), KC.collision_distribution_coefficient(collision))
        for (term, _) in expression
            term in seen && continue
            push!(seen, term)
            KC.finite_width_reduction_kind(term) === KC.FiniteWidthConvolution || continue
            polynomial = skeleton_external_response_polynomial(collision, term)
            iszero(KC.evaluate_occupation_polynomial(polynomial, background)) && continue
            linearization = KC.evaluate_occupation_linearization(
                KC.occupation_linearization(polynomial), background
            )
            iszero(linearization) && continue
            return term, polynomial, linearization
        end
    end
    return error("skeleton collision contains no usable convolution response term")
end

function skeleton_external_response_lines(collision)
    lines = KC.SpectralLineIdentity{Boson}[]
    target = KC.target_family(collision)
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
            else
                continue
            end
            push!(lines, KC._finite_width_external_spectral_line(term, target))
        end
    end
    unique!(lines)
    return lines
end

function skeleton_external_scalar_linewidth(kernel, occupation_value, γ_value)
    background = KC.OccupationBackground(_ -> occupation_value, Rational{Int})
    evaluated = KC.evaluate_spectral_self_energy_background(kernel, background)
    sector, coefficient = only(KC.background_spectral_self_energy_terms(evaluated))

    # Benchmark-local closure: the package deliberately leaves the routed loop integral explicit.
    # With ∫_q 1 = 1, the generated first-order source gives Γ = 4γ n and δΓ = 4γ δn.
    unit_normalized_loop_integral = (_, value) -> real(value)
    return γ_value * unit_normalized_loop_integral(sector, coefficient)
end

@testset "generated skeleton loss closes the full external finite-width response" begin
    ordinary_collision = skeleton_external_response_collision(; skeleton=false)
    skeleton_collision = skeleton_external_response_collision(; skeleton=true)

    ordinary_kinds = skeleton_external_response_kinds(ordinary_collision)
    skeleton_kinds = skeleton_external_response_kinds(skeleton_collision)
    @test KC.FiniteWidthFactorized in ordinary_kinds
    @test KC.FiniteWidthConvolution in ordinary_kinds
    @test !(KC.FiniteWidthFactorized in skeleton_kinds)
    @test KC.FiniteWidthConvolution in skeleton_kinds

    linewidth_kernel = skeleton_external_response_linewidth_kernel()
    γ_value = 1 // 4
    background_value = 3 // 1
    unit_variation = 1 // 1
    Γbar = skeleton_external_scalar_linewidth(linewidth_kernel, background_value, γ_value)
    δΓ = skeleton_external_scalar_linewidth(linewidth_kernel, unit_variation, γ_value)
    @test Γbar == 3 // 1
    @test δΓ == 1 // 1

    background = KC.OccupationBackground(_ -> background_value, Rational{Int})
    term, polynomial, background_linearization = skeleton_external_response_convolution(
        skeleton_collision, background
    )
    reduction = KC.finite_width_convolution_reduction(term)
    target = KC.target_family(skeleton_collision)
    external_line = KC._finite_width_external_spectral_line(term, target)
    internal_lines = vcat(
        KC.convolution_pivot_lines(reduction), [KC.convolution_dependent_line(reduction)]
    )
    @test length(internal_lines) == 3
    @test length(unique(vcat(internal_lines, [external_line]))) == 4
    @test abs(KC.convolution_external_coefficient(reduction)) == 1 // 1
    @test all(
        coefficient -> abs(coefficient) == 1 // 1, KC.convolution_coefficients(reduction)
    )

    lines = skeleton_external_response_lines(skeleton_collision)
    model = KC.LorentzianSpectralModel(
        Dict(line => KC.LorentzianSpectralData(0 // 1, Γbar, 1 // 1) for line in lines)
    )
    model_variation = KC.LorentzianSpectralModelVariation(
        Dict(
            line => KC.LorentzianSpectralDataVariation(0 // 1, δΓ, 0 // 1) for line in lines
        ),
    )
    atom, δP = first(KC.background_occupation_linearization_terms(background_linearization))
    Pbar = KC.evaluate_occupation_polynomial(polynomial, background)
    response = KC.OccupationSpectralResponse(Dict(atom => model_variation))

    full = @inferred KC.linearize_external_spectral_collision(
        skeleton_collision, model, background, response
    )

    J = KC.convolution_jacobian(reduction)
    expected_weight = J / Γbar
    expected_δweight = -J * δΓ / Γbar^2
    @test KC.external_spectral_projection_linewidth(reduction, model, external_line) ==
        4 * Γbar
    @test KC.external_spectral_projection_mismatch(reduction, model, external_line) ==
        0 // 1
    @test KC.evaluate_external_spectral_projection(reduction, model, external_line) ==
        expected_weight
    @test KC.evaluate_external_spectral_projection_variation(
        reduction, model, model_variation, external_line
    ) == expected_δweight

    occupation_terms = Dict(
        KC.background_occupation_linearization_terms(
            KC.external_spectral_occupation_response_terms(full)[term]
        ),
    )
    spectral_terms = Dict(
        KC.background_occupation_linearization_terms(
            KC.external_spectral_model_response_terms(full)[term]
        ),
    )
    full_terms = Dict(
        KC.background_occupation_linearization_terms(
            KC.external_spectral_linearized_terms(full)[term]
        ),
    )

    @test occupation_terms[atom] == expected_weight * δP
    @test spectral_terms[atom] == Pbar * expected_δweight
    @test full_terms[atom] == occupation_terms[atom] + spectral_terms[atom]
    @test full_terms[atom] != occupation_terms[atom]
end
