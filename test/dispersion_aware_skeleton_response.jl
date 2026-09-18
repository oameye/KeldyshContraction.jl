using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields dispersion_aware_skeleton_ϕ::Boson

function dispersion_aware_loss_lagrangian()
    c = dispersion_aware_skeleton_ϕ[Classical]
    q = dispersion_aware_skeleton_ϕ[Quantum]
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

function dispersion_aware_kinetic(order::Val{O}, edges::Val{E}) where {O,E}
    L = dispersion_aware_loss_lagrangian()
    G = DressedPropagator(L, order, edges; simplify=true, _set_reg_to_zero=false)
    ΣF = KC.skeleton_self_energy(fourier_transform(G))
    return kinetic_expression(wigner_transform(ΣF; gradient_order=Val(0)))
end

function dispersion_aware_collision()
    kinetic = dispersion_aware_kinetic(Val(2), Val(5))
    offshell = off_shell_collision_expression(kinetic)
    spectral = spectral_dispersive_collision(offshell)
    return KC.finite_width_frequency_collision(spectral)
end

function dispersion_aware_linewidth(occupation_value, γ_value)
    kernel = KC.spectral_self_energy_kernel(dispersion_aware_kinetic(Val(1), Val(3)))
    background = KC.OccupationBackground(_ -> occupation_value, Rational{Int})
    evaluated = KC.evaluate_spectral_self_energy_background(kernel, background)
    _, coefficient = only(KC.background_spectral_self_energy_terms(evaluated))
    return γ_value * real(coefficient)
end

function dispersion_aware_convolution(collision)
    background = KC.OccupationBackground(_ -> 1 // 1, Rational{Int})
    target = KC.target_family(collision)
    seen = Set{Any}()
    for (expression, include_external) in (
        KC.collision_offset(collision) => false,
        KC.collision_distribution_coefficient(collision) => true,
    )
        for (term, coefficient) in expression
            term in seen && continue
            push!(seen, term)
            KC.finite_width_reduction_kind(term) === KC.FiniteWidthConvolution || continue
            polynomial = KC._finite_width_occupation_polynomial(
                term, coefficient, target, include_external
            )
            iszero(KC.evaluate_occupation_polynomial(polynomial, background)) && continue
            return term
        end
    end
    return error("generated skeleton collision contains no usable convolution term")
end

function physical_momentum(momentum, values)
    x = 0 // 1
    y = 0 // 1
    for i in eachindex(values)
        coefficient = momentum[i]
        x += coefficient * values[i][1]
        y += coefficient * values[i][2]
    end
    return x, y
end

function quadratic_energy(line, values)
    px, py = physical_momentum(KC.spectral_line_momentum(line), values)
    return (px^2 + py^2) / 2
end

@testset "dispersion-aware four-line skeleton spectral geometry" begin
    collision = dispersion_aware_collision()
    term = dispersion_aware_convolution(collision)
    reduction = KC.finite_width_convolution_reduction(term)
    target = KC.target_family(collision)
    external_line = KC._finite_width_external_spectral_line(term, target)
    internal_lines = vcat(
        KC.convolution_pivot_lines(reduction), [KC.convolution_dependent_line(reduction)]
    )
    lines = unique(vcat(internal_lines, [external_line]))
    @test length(lines) == 4

    source = KC.finite_width_source_term(term)
    basis = KC.momentum_basis(source)
    external = KC.external_wigner_momentum(source)
    external_index = only(i for (i, variable) in enumerate(basis) if variable == external)
    loop_indices = [i for i in eachindex(basis.variables) if i != external_index]
    @test length(loop_indices) == 2

    values = [(0 // 1, 0 // 1) for _ in eachindex(basis.variables)]
    values[external_index] = (1 // 1, 0 // 1)
    values[loop_indices[1]] = (0 // 1, 1 // 1)
    values[loop_indices[2]] = (1 // 1, 1 // 1)

    energies = sort([quadratic_energy(line, values) for line in lines])
    @test energies == [1 // 2, 1 // 2, 1 // 1, 2 // 1]

    γ_value = 1 // 16
    Γbar = dispersion_aware_linewidth(1 // 1, γ_value)
    δΓ = dispersion_aware_linewidth(1 // 1, γ_value)
    @test Γbar == 1 // 4
    @test δΓ == 1 // 4

    model = KC.LorentzianSpectralModel(
        Dict(
            line => KC.LorentzianSpectralData(quadratic_energy(line, values), Γbar, 1 // 1) for
            line in lines
        ),
    )
    variation = KC.LorentzianSpectralModelVariation(
        Dict(
            line => KC.LorentzianSpectralDataVariation(0 // 1, δΓ, 0 // 1) for
            line in lines
        ),
    )

    mismatch = KC.external_spectral_projection_mismatch(reduction, model, external_line)
    width = KC.external_spectral_projection_linewidth(reduction, model, external_line)
    J = KC.convolution_jacobian(reduction)

    @test abs(mismatch) == 1 // 1
    @test width == 1 // 1
    @test KC.evaluate_external_spectral_projection(reduction, model, external_line) == 4J / 5
    @test KC.evaluate_external_spectral_projection_variation(
        reduction, model, variation, external_line
    ) == 12J / 25

    # Off resonance, increasing a narrow linewidth initially increases spectral overlap.
    # This is the opposite sign from the resonant equal-energy compiler oracle and is why
    # the physical Maxwell closure must retain the kinetic-energy mismatch pointwise.
    @test KC.evaluate_external_spectral_projection_variation(
        reduction, model, variation, external_line
    ) > 0
end
