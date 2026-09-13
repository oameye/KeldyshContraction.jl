using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields generated_kernel_ϕ::Boson

function generated_kernel_elastic_collision()
    c = generated_kernel_ϕ[Classical]
    q = generated_kernel_ϕ[Quantum]
    interaction =
        -(1 // 2) * ((c^2 + q^2) * bar(c) * bar(q) + c * q * (bar(c)^2 + bar(q)^2))
    L = InteractionLagrangian(interaction, :g)
    G = DressedPropagator(L, Val(2), Val(5); simplify=true)
    GF = fourier_transform(G)
    ΣW = wigner_transform(SelfEnergy(GF); gradient_order=Val(0))
    return spectral_dispersive_collision(
        off_shell_collision_expression(kinetic_expression(ΣW))
    )
end

function generated_kernel_loss_collision(order::Val{O}, edges::Val{E}) where {O,E}
    c = generated_kernel_ϕ[Classical]
    q = generated_kernel_ϕ[Quantum]
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
    G = DressedPropagator(L, order, edges; simplify=true, _set_reg_to_zero=false)
    GF = fourier_transform(G)
    ΣW = wigner_transform(SelfEnergy(GF); gradient_order=Val(0))
    return spectral_dispersive_collision(
        off_shell_collision_expression(kinetic_expression(ΣW))
    )
end

function generated_kernel_mixed_collision()
    c = generated_kernel_ϕ[Classical]
    q = generated_kernel_ϕ[Quantum]
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus

    elastic = -(1 // 2) * ((c^2 + q^2) * bar(c) * bar(q) + c * q * (bar(c)^2 + bar(q)^2))
    loss =
        (1 // 2) *
        im *
        (
            bar(c) * bar(q) * (c(minus)^2 + q(minus)^2) -
            c(plus) * q(plus) * (bar(c)^2 + bar(q)^2) +
            2 * bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
        )

    L = InteractionLagrangian(elastic, :g) + InteractionLagrangian(loss, :γ)
    G = DressedPropagator(L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false)
    mixed_parameter = KC.ParameterMonomial(:g) * KC.ParameterMonomial(:γ)
    GF = fourier_transform(G[mixed_parameter])
    ΣW = wigner_transform(SelfEnergy(GF); gradient_order=Val(0))
    return spectral_dispersive_collision(
        off_shell_collision_expression(kinetic_expression(ΣW))
    )
end

function generated_regular_kernel(collision)
    reduced = @inferred reduce_frequency_collision(collision)
    occupation = @inferred occupation_reduced_expression(reduced)
    quotient = @inferred quotient_loop_momenta(occupation)
    kernel = @inferred collision_kernel(quotient)
    return reduced, occupation, quotient, kernel
end

function generated_expected_kernel(occupation, sector, polynomial)
    expression = typeof(occupation)(
        Dict(sector => polynomial),
        target_family(occupation),
        parameters(occupation),
        wigner_context(occupation),
    )
    return @inferred collision_kernel(expression)
end

function generated_kernel_basis(sector)
    basis = KC.momentum_basis(sector)
    external = external_wigner_momentum(sector)
    external_index = only(i for (i, variable) in enumerate(basis) if variable == external)
    loop_indices = Int[i for i in eachindex(basis.variables) if i != external_index]
    return basis, external_index, loop_indices
end

@testset "generated first-order loss CollisionKernel" begin
    collision = generated_kernel_loss_collision(Val(1), Val(3))
    reduced, occupation, quotient, kernel = generated_regular_kernel(collision)

    @test isempty(reduced_dependent_terms(reduced))
    @test isempty(reduced_causal_terms(reduced))
    @test isempty(reduced_trotter_terms(reduced))
    @test length(occupation_reduced_terms(occupation)) == 1
    @test length(loop_quotient_terms(quotient)) == 1

    sector, _ = only(occupation_reduced_terms(occupation))
    basis, external_index, loop_indices = generated_kernel_basis(sector)
    @test length(loop_indices) == 1
    k = KC.basis_momentum(basis, external_index)
    q = KC.basis_momentum(basis, only(loop_indices))
    C = KC.ComplexRationals
    n(momentum) = OccupationPolynomial(OccupationAtom(generated_kernel_ϕ, momentum), one(C))
    expected = -4 * n(k) * n(q)

    expected_kernel = generated_expected_kernel(occupation, sector, expected)
    @test collision_kernel_terms(kernel) == collision_kernel_terms(expected_kernel)
    @test parameters(kernel) == KC.ParameterMonomial(:γ)
end

@testset "generated elastic g² quotient equals standard Bose kernel" begin
    collision = generated_kernel_elastic_collision()
    reduced, occupation, _, kernel = generated_regular_kernel(collision)

    @test isempty(reduced_dependent_terms(reduced))
    @test isempty(reduced_causal_terms(reduced))
    @test isempty(reduced_trotter_terms(reduced))
    sector, _ = only(occupation_reduced_terms(occupation))
    @test parameters(sector) == KC.ParameterMonomial(:g)^2

    basis, external_index, loop_indices = generated_kernel_basis(sector)
    @test length(loop_indices) == 2
    k = KC.basis_momentum(basis, external_index)
    q = KC.basis_momentum(basis, loop_indices[1])
    r = KC.basis_momentum(basis, loop_indices[2])
    p = -k + q + r
    C = KC.ComplexRationals
    n(momentum) = OccupationPolynomial(OccupationAtom(generated_kernel_ϕ, momentum), one(C))
    nk, np, nq, nr = n(k), n(p), n(q), n(r)
    one_n = one(typeof(nk))
    standard =
        2 * ((one_n + nk) * (one_n + np) * nq * nr - nk * np * (one_n + nq) * (one_n + nr))

    expected_kernel = generated_expected_kernel(occupation, sector, standard)
    @test collision_kernel_terms(kernel) == collision_kernel_terms(expected_kernel)
    @test all(
        key ->
            length(frequency_support(key).shells) == 1 &&
            isempty(frequency_support(key).principal_values),
        keys(collision_kernel_terms(kernel)),
    )
end

@testset "generated mixed gγ PV quotient" begin
    collision = generated_kernel_mixed_collision()
    reduced, occupation, _, kernel = generated_regular_kernel(collision)

    @test isempty(reduced_dependent_terms(reduced))
    @test isempty(reduced_causal_terms(reduced))
    @test isempty(reduced_trotter_terms(reduced))
    sector, _ = only(occupation_reduced_terms(occupation))
    @test parameters(sector) == KC.ParameterMonomial(:g) * KC.ParameterMonomial(:γ)
    @test isempty(frequency_support(sector).shells)
    @test length(frequency_support(sector).principal_values) == 1

    basis, external_index, loop_indices = generated_kernel_basis(sector)
    @test length(loop_indices) == 2
    k = KC.basis_momentum(basis, external_index)
    q = KC.basis_momentum(basis, loop_indices[1])
    r = KC.basis_momentum(basis, loop_indices[2])
    p = -k + q + r
    C = KC.ComplexRationals
    n(momentum) = OccupationPolynomial(OccupationAtom(generated_kernel_ϕ, momentum), one(C))
    frozen_309 = -8 * n(k) * n(p) - 16 * n(k) * n(p) * n(q)

    expected_kernel = generated_expected_kernel(occupation, sector, frozen_309)
    @test collision_kernel_terms(kernel) == collision_kernel_terms(expected_kernel)
    @test all(
        key ->
            isempty(frequency_support(key).shells) &&
            length(frequency_support(key).principal_values) == 1,
        keys(collision_kernel_terms(kernel)),
    )
end

@testset "generated γ² regular quotient equals independent 2𝒩 kernel" begin
    collision = generated_kernel_loss_collision(Val(2), Val(5))
    reduced, occupation, _, kernel = generated_regular_kernel(collision)

    @test !isempty(reduced_dependent_terms(reduced)) ||
        !isempty(reduced_causal_terms(reduced))
    sector, _ = only(occupation_reduced_terms(occupation))
    @test parameters(sector) == KC.ParameterMonomial(:γ)^2

    basis, external_index, loop_indices = generated_kernel_basis(sector)
    @test length(loop_indices) == 2
    k = KC.basis_momentum(basis, external_index)
    q1 = KC.basis_momentum(basis, loop_indices[1])
    q2 = KC.basis_momentum(basis, loop_indices[2])
    q = -k + q1 + q2
    C = KC.ComplexRationals
    n(momentum) = OccupationPolynomial(OccupationAtom(generated_kernel_ϕ, momentum), one(C))
    nk, nq, nq1, nq2 = n(k), n(q), n(q1), n(q2)
    standard =
        2 * (
            nq1 * nq2 +
            nq1 * nq2 * nq +
            nq1 * nq2 * nk +
            nq * nk +
            nq1 * nq * nk +
            nq2 * nq * nk
        )

    expected_kernel = generated_expected_kernel(occupation, sector, standard)
    @test collision_kernel_terms(kernel) == collision_kernel_terms(expected_kernel)
    @test all(
        key ->
            length(frequency_support(key).shells) == 1 &&
            isempty(frequency_support(key).principal_values),
        keys(collision_kernel_terms(kernel)),
    )
end
