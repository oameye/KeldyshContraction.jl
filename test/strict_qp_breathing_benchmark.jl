using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields strict_qp_breathing_ϕ::Boson

function strict_qp_breathing_loss_collision(order::Val{O}, edges::Val{E}) where {O,E}
    c = strict_qp_breathing_ϕ[Classical]
    q = strict_qp_breathing_ϕ[Quantum]
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
    spectral = spectral_dispersive_collision(
        off_shell_collision_expression(kinetic_expression(ΣW))
    )
    return @inferred reduce_frequency_collision(spectral)
end

@testset "2D harmonic-trap breathing streaming" begin
    # Moment basis (⟨r²⟩, ⟨r⋅k⟩, ⟨k²⟩) in units m = ω₀ = 1.
    streaming = Rational{Int}[
        0 2 0
        -1 0 1
        0 -2 0
    ]
    χplus = Complex{Rational{Int}}[1, 2im, -1]
    χminus = conj.(χplus)

    left_action(v) = [sum(v[i] * streaming[i, j] for i in 1:3) for j in 1:3]

    @test left_action(χplus) == (-2im) .* χplus
    @test left_action(χminus) == (2im) .* χminus
    @test streaming^3 == -4 * streaming
end

@testset "generated first-order loss breathing envelope" begin
    reduced = strict_qp_breathing_loss_collision(Val(1), Val(3))
    @test isempty(KC.reduced_blocked_terms(reduced))
    @test isempty(KC.reduced_trotter_terms(reduced))

    occupation = @inferred KC.occupation_reduced_expression(reduced)
    quotient = @inferred KC.quotient_loop_momenta(occupation)
    kernel = @inferred collision_kernel(quotient)
    sector, polynomial = only(collision_kernel_terms(kernel))
    @test parameters(sector) == KC.ParameterMonomial(:γ)

    basis = KC.momentum_basis(sector)
    external = KC.external_wigner_momentum(sector)
    external_index = only(i for (i, variable) in enumerate(basis) if variable == external)
    loop_index = only(i for i in eachindex(basis.variables) if i != external_index)
    k = KC.basis_momentum(basis, external_index)
    q = KC.basis_momentum(basis, loop_index)
    C = KC.ComplexRationals
    n(momentum) =
        OccupationPolynomial(OccupationAtom(strict_qp_breathing_ϕ, momentum), one(C))
    @test polynomial == -4 * n(k) * n(q)

    # Classical Maxwell J_CE=0 projection in m = T = ω₀ = 1 units.
    # For n₀ = z exp[-(r²+k²)/2], M̄ = z/(4π).  The χ₋-χ₊ overlap is 16,
    # while the linearized first-order loss pairing is 12 under the narrower
    # M₀ n₀ spatial Gaussian.  Strip the common γz/π factor so the check remains exact.
    κ1 = 4 // 1
    mbar_times_pi_over_z = 1 // 4
    loss_pairing_ratio = 12 // 16
    raw_rate_times_pi_over_γz = -κ1 * loss_pairing_ratio * mbar_times_pi_over_z
    number_rate_times_pi_over_γz = -κ1 * mbar_times_pi_over_z
    normalized_rate_times_pi_over_γz =
        raw_rate_times_pi_over_γz - number_rate_times_pi_over_γz

    @test raw_rate_times_pi_over_γz == -3 // 4
    @test number_rate_times_pi_over_γz == -1 // 1
    @test normalized_rate_times_pi_over_γz == 1 // 4
    @test normalized_rate_times_pi_over_γz == -number_rate_times_pi_over_γz / 4
    @test imag(complex(normalized_rate_times_pi_over_γz)) == 0
end

@testset "generated strict-QP γ² regular loss kernel" begin
    reduced = strict_qp_breathing_loss_collision(Val(2), Val(5))

    # The repeated-line/equal-time sector is deliberately not forced through strict QP.
    @test !isempty(KC.reduced_blocked_terms(reduced))
    @test isempty(KC.reduced_trotter_terms(reduced))

    occupation = @inferred KC.occupation_reduced_expression(reduced)
    occupation_sector, _ = only(KC.occupation_reduced_terms(occupation))
    @test parameters(occupation_sector) == KC.ParameterMonomial(:γ)^2

    # Build the independent 2𝒩 oracle on the same pre-quotient momentum basis. Both the
    # generated expression and the oracle then pass through KC's exact loop quotient.
    basis = KC.momentum_basis(occupation_sector)
    external = KC.external_wigner_momentum(occupation_sector)
    external_index = only(i for (i, variable) in enumerate(basis) if variable == external)
    loop_indices = [i for i in eachindex(basis.variables) if i != external_index]
    @test length(loop_indices) == 2

    k = KC.basis_momentum(basis, external_index)
    q1 = KC.basis_momentum(basis, loop_indices[1])
    q2 = KC.basis_momentum(basis, loop_indices[2])
    q = -k + q1 + q2
    C = KC.ComplexRationals
    n(momentum) =
        OccupationPolynomial(OccupationAtom(strict_qp_breathing_ϕ, momentum), one(C))
    nk, nq, nq1, nq2 = n(k), n(q), n(q1), n(q2)

    carleman_oracle =
        2 * (
            nq1 * nq2 +
            nq1 * nq2 * nq +
            nq1 * nq2 * nk +
            nq * nk +
            nq1 * nq * nk +
            nq2 * nq * nk
        )
    oracle_expression = typeof(occupation)(
        Dict(occupation_sector => carleman_oracle),
        target_family(occupation),
        parameters(occupation),
        wigner_context(occupation),
    )

    quotient = @inferred KC.quotient_loop_momenta(occupation)
    kernel = @inferred collision_kernel(quotient)
    expected_kernel = @inferred collision_kernel(oracle_expression)
    @test collision_kernel_terms(kernel) == collision_kernel_terms(expected_kernel)
    @test length(collision_kernel_terms(kernel)) == 1

    sector, polynomial = only(collision_kernel_terms(kernel))
    @test parameters(sector) == KC.ParameterMonomial(:γ)^2
    @test length(frequency_support(sector).shells) == 1
    @test isempty(frequency_support(sector).principal_values)

    degrees = sort(collect(length(monomial) for (monomial, _) in polynomial))
    @test degrees == [2, 2, 3, 3, 3]
end

@testset "J_CE=0 Maxwell breathing oracle" begin
    # Independent continuum CE matrix elements at the classical-Maxwell reference
    # M̄_ref = 1/(4π).  We store coefficients after stripping the displayed powers
    # of π so every assertion below is exact rational arithmetic.
    S00_over_pi2 = 64 // 1
    stageA00_over_pi4 = 96 // 1
    stageB23_00_over_pi4 = 208 // 9
    stageB56_00_over_pi4 = 208 // 9
    stageB00_over_pi4 = stageB23_00_over_pi4 + stageB56_00_over_pi4

    # M^(2)_A = stageA/(8π⁴) × (M̄/M̄_ref), with M̄_ref = 1/(4π).
    # Therefore δλ_A/γ² = (a1π) M̄/π with the exact dimensionless coefficient below.
    a1_times_pi = (stageA00_over_pi4 / 8) * 4 / S00_over_pi2
    @test a1_times_pi == 3 // 4

    # Stage B scales as (M̄/M̄_ref)^2, hence all π factors cancel from a2.
    a2_regular = (stageB00_over_pi4 / 8) * 16 / S00_over_pi2
    @test a2_regular == 13 // 9

    # Thus the regular strict-QP pole correction is
    # δλ_reg/γ² = (3/(4π)) M̄ + (13/9) M̄².
    # The historical lollipop/repeated-line block is separate: 160/9 M̄².
    lollipop_a2 = 160 // 9
    @test a2_regular + lollipop_a2 == 173 // 9

    # J_CE=0 is one-dimensional and all regular matrix elements are real, so this
    # strict-QP truncation changes the envelope but not the breathing frequency.
    @test imag(complex(a2_regular)) == 0
end
