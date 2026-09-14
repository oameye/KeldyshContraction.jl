using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields occupation_loss_ϕ::Boson

function generated_loss_occupation_collision()
    c = occupation_loss_ϕ[Classical]
    q = occupation_loss_ϕ[Quantum]
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
    GF = fourier_transform(G)
    ΣW = wigner_transform(SelfEnergy(GF); gradient_order=Val(0))
    spectral = spectral_dispersive_collision(
        off_shell_collision_expression(kinetic_expression(ΣW))
    )
    return KC.canonical_frequency_collision(spectral)
end

@testset "generated γ² raw occupation reduction" begin
    canonical = generated_loss_occupation_collision()
    reduced = @inferred KC.reduce_frequency_collision(canonical)

    @test isempty(KC.reduced_trotter_terms(reduced))
    @test !isempty(KC.reduced_blocked_terms(reduced))
    @test length(KC.reduced_regular_terms(reduced)) == 1

    occupation = @inferred KC.occupation_reduced_expression(reduced)
    @test length(KC.occupation_reduced_terms(occupation)) == 1
    sector, polynomial = only(KC.occupation_reduced_terms(occupation))
    @test parameters(sector) == KC.ParameterMonomial(:γ)^2
    @test length(KC.frequency_support(sector).shells) == 1
    @test isempty(KC.frequency_support(sector).principal_values)

    basis = KC.momentum_basis(sector)
    external_variable = KC.external_wigner_momentum(sector)
    external_index = only(
        i for (i, variable) in enumerate(basis) if variable == external_variable
    )
    loop_indices = [i for i in eachindex(basis.variables) if i != external_index]
    @test length(loop_indices) == 2

    k = KC.basis_momentum(basis, external_index)
    q1 = KC.basis_momentum(basis, loop_indices[1])
    q2 = KC.basis_momentum(basis, loop_indices[2])
    q = -k + q1 + q2

    coefficient = one(KC.ComplexRationals)
    n(momentum) =
        KC.OccupationPolynomial(KC.OccupationAtom(occupation_loss_ϕ, momentum), coefficient)
    nk, nq, nq1, nq2 = n(k), n(q), n(q1), n(q2)

    symmetric_oracle =
        2 * (
            nq1 * nq2 +
            nq1 * nq2 * nq +
            nq1 * nq2 * nk +
            nq * nk +
            nq1 * nq * nk +
            nq2 * nq * nk
        )

    raw_oracle =
        nq * nq2 + 2 * nq1 * nq2 * nq - nq1 * nq +
        4 * nq1 * nq * nk +
        2 * nq * nk +
        nq2 +
        2 * nq1 * nq2 +
        2 * nq1 * nq2 * nk +
        nq2 * nk - nq1 - nq1 * nk
    swapped_raw_oracle =
        nq * nq1 + 2 * nq1 * nq2 * nq - nq2 * nq +
        4 * nq2 * nq * nk +
        2 * nq * nk +
        nq1 +
        2 * nq1 * nq2 +
        2 * nq1 * nq2 * nk +
        nq1 * nk - nq2 - nq2 * nk

    @test polynomial == raw_oracle
    @test polynomial != symmetric_oracle
    @test polynomial + swapped_raw_oracle == 2 * symmetric_oracle
end
