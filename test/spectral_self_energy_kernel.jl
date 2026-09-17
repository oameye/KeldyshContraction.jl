using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields spectral_width_ϕ::Boson

function generated_first_order_loss_self_energy()
    c = spectral_width_ϕ[Classical]
    q = spectral_width_ϕ[Quantum]
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
    ΣF = SelfEnergy(fourier_transform(G))
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    return @inferred kinetic_expression(ΣW)
end

@testset "generated first-order loss spectral self-energy" begin
    Σ = generated_first_order_loss_self_energy()
    kernel = @inferred KC.spectral_self_energy_kernel(Σ)

    @test KC.order(kernel) == 1
    @test KC.statistics(kernel) === Boson
    @test KC.parameters(kernel) == KC.ParameterMonomial(:γ)
    @test KC.target_family(kernel) === spectral_width_ϕ
    @test KC.gradient_order(kernel) == Val(0)
    @test isempty(KC.spectral_self_energy_blocked_terms(kernel))
    @test isempty(KC.spectral_self_energy_trotter_terms(kernel))
    @test length(KC.spectral_self_energy_terms(kernel)) == 1

    sector, polynomial = only(KC.spectral_self_energy_terms(kernel))
    @test KC.parameters(sector) == KC.ParameterMonomial(:γ)
    @test isempty(KC.frequency_support(sector).shells)
    @test isempty(KC.frequency_support(sector).principal_values)

    basis = KC.momentum_basis(sector)
    external = KC.external_wigner_momentum(sector)
    external_index = only(i for (i, variable) in enumerate(basis) if variable == external)
    loop_index = only(i for i in eachindex(basis.variables) if i != external_index)
    q = KC.basis_momentum(basis, loop_index)
    C = KC.ComplexRationals
    nq = KC.OccupationPolynomial(KC.OccupationAtom(spectral_width_ϕ, q), one(C))

    # A_Σ = im(Σᴿ-Σᴬ) = Γ.  Equal-time regularisation cancels the vacuum piece,
    # leaving Γ(k) = 4γ ∫_q n_q for bosonic two-body loss.
    @test polynomial == 4 * nq
end

@testset "spectral self-energy evaluates occupation background before momentum integration" begin
    kernel = @inferred KC.spectral_self_energy_kernel(
        generated_first_order_loss_self_energy()
    )
    source_sector, _ = only(KC.spectral_self_energy_terms(kernel))
    background = KC.OccupationBackground(_ -> 3 // 1, Rational{Int64})

    evaluated = @inferred KC.evaluate_spectral_self_energy_background(kernel, background)
    @test evaluated isa KC.BackgroundSpectralSelfEnergyKernel
    @test KC.order(evaluated) == KC.order(kernel)
    @test KC.statistics(evaluated) === KC.statistics(kernel)
    @test KC.target_family(evaluated) === KC.target_family(kernel)
    @test KC.parameters(evaluated) == KC.parameters(kernel)
    @test KC.gradient_order(evaluated) == KC.gradient_order(kernel)
    @test KC.wigner_context(evaluated) == KC.wigner_context(kernel)
    @test isempty(KC.spectral_self_energy_blocked_terms(evaluated))
    @test isempty(KC.spectral_self_energy_trotter_terms(evaluated))

    sector, coefficient = only(KC.background_spectral_self_energy_terms(evaluated))
    @test sector == source_sector
    @test coefficient == 12 // 1

    # This is the coefficient of the still-explicit q integral.  No hidden momentum
    # integration is performed by background evaluation.
    basis = KC.momentum_basis(sector)
    external = KC.external_wigner_momentum(sector)
    @test count(variable -> variable != external, basis) == 1

    # A downstream calculation supplies its own momentum integral.  With a benchmark-local
    # normalized rule ∫_q 1 = 1, the resulting scalar Γ feeds the existing Lorentzian API
    # directly; KC does not need another integration-wrapper abstraction.
    unit_normalized_loop_integral = (_, value) -> real(value)
    linewidth = unit_normalized_loop_integral(sector, coefficient)
    data = KC.LorentzianSpectralData(2 // 1, linewidth)
    @test KC.spectral_linewidth(data) == 12 // 1
    @test KC.spectral_squared_weight(data) == 1 // 6

    zero_background = KC.OccupationBackground(_ -> 0 // 1, Rational{Int64})
    zero_evaluated = @inferred KC.evaluate_spectral_self_energy_background(
        kernel, zero_background
    )
    @test isempty(KC.background_spectral_self_energy_terms(zero_evaluated))
end
