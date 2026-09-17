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
