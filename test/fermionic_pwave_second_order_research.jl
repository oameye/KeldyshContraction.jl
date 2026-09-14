using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields second_order_pwave_ψ::Fermion

function pwave_branch_pairs(; shifted::Bool)
    ψ₁ = second_order_pwave_ψ[One]
    ψ₂ = second_order_pwave_ψ[Two]
    bψ₁ = bar(ψ₁)
    bψ₂ = bar(ψ₂)

    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus

    ψplus = shifted ? ψ₁(minus) + ψ₂(minus) : ψ₁ + ψ₂
    ∂ψplus = if shifted
        partial(ψ₁(minus), :x) + partial(ψ₂(minus), :x)
    else
        partial(ψ₁, :x) + partial(ψ₂, :x)
    end
    ψplus_forward = shifted ? ψ₁(plus) + ψ₂(plus) : ψ₁ + ψ₂
    ∂ψplus_forward = if shifted
        partial(ψ₁(plus), :x) + partial(ψ₂(plus), :x)
    else
        partial(ψ₁, :x) + partial(ψ₂, :x)
    end
    ψminus = shifted ? ψ₁(plus) - ψ₂(plus) : ψ₁ - ψ₂
    ∂ψminus = if shifted
        partial(ψ₁(plus), :x) - partial(ψ₂(plus), :x)
    else
        partial(ψ₁, :x) - partial(ψ₂, :x)
    end

    bψplus = bψ₁ + bψ₂
    bψminus = bψ₂ - bψ₁
    ∂bψplus = partial(bψ₁, :x) + partial(bψ₂, :x)
    ∂bψminus = partial(bψ₂, :x) - partial(bψ₁, :x)

    Pplus = ψplus * ∂ψplus
    Pplus_forward = ψplus_forward * ∂ψplus_forward
    Pminus = ψminus * ∂ψminus
    Pplus_dagger = ∂bψplus * bψplus
    Pminus_dagger = ∂bψminus * bψminus
    return Pplus, Pplus_forward, Pminus, Pplus_dagger, Pminus_dagger
end

function fermionic_pwave_elastic_lagrangian()
    Pplus, _, Pminus, Pplus_dagger, Pminus_dagger = pwave_branch_pairs(; shifted=false)

    # H_int = (g_p/2) P†P with P = ψ∂xψ. The LO branch sums each carry the
    # 1/√2 rotation factors, leaving the exact branch coefficient -1/8.
    elastic = -(1 // 8) * (Pplus_dagger * Pplus - Pminus_dagger * Pminus)
    return InteractionLagrangian(elastic, :gp)
end

function fermionic_pwave_loss_lagrangian_second_order()
    Pplus_minus, Pplus_plus, Pminus_plus, Pplus_dagger, Pminus_dagger = pwave_branch_pairs(;
        shifted=true
    )

    # Same finite-Trotter Lindblad convention as the certified O(γp) oracle in #312.
    loss =
        (1 // 8) *
        im *
        (
            (Pplus_dagger - Pminus_dagger) * Pplus_minus -
            (Pplus_plus - Pminus_plus) * Pminus_dagger
        )
    return InteractionLagrangian(loss, :γp)
end

function generated_second_order_pipeline(G, parameter)
    GF = fourier_transform(G[parameter])
    ΣF = SelfEnergy(GF)
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    kinetic = kinetic_expression(ΣW)
    off_shell = off_shell_collision_expression(kinetic)
    spectral = spectral_dispersive_collision(off_shell)
    canonical = KC.canonical_frequency_collision(spectral)
    reduced = KC.reduce_frequency_collision(canonical)
    occupation = KC.occupation_reduced_expression(reduced)
    quotient = KC.quotient_loop_momenta(occupation)
    kernel = KC.collision_kernel(quotient)
    return (;
        GF,
        ΣF,
        ΣW,
        kinetic,
        off_shell,
        spectral,
        canonical,
        reduced,
        occupation,
        quotient,
        kernel,
    )
end

function research_second_order_census(result)
    spectral_terms =
        sum(1 for _ in collision_offset(result.spectral)) +
        sum(1 for _ in collision_distribution_coefficient(result.spectral))
    canonical_regular = length(KC.canonical_frequency_expressions(result.canonical))
    canonical_shifted = length(KC.shifted_frequency_terms(result.canonical))
    regular = length(KC.reduced_regular_terms(result.reduced))
    blocked = length(KC.reduced_blocked_terms(result.reduced))
    trotter = length(KC.reduced_trotter_terms(result.reduced))
    occupation = length(KC.occupation_reduced_terms(result.occupation))
    quotient = length(KC.loop_quotient_terms(result.quotient))
    kernel = length(KC.collision_kernel_terms(result.kernel))
    return (;
        spectral_terms,
        canonical_regular,
        canonical_shifted,
        regular,
        blocked,
        trotter,
        occupation,
        quotient,
        kernel,
    )
end

function pwave_axis_polynomial(momentum)
    C = KC.ComplexRationals
    component = KC.MomentumComponent(momentum, :x)
    monomial = KC.MomentumMonomial(KC.MomentumComponent[component])
    return KC.MomentumPolynomial(monomial, one(C))
end

function certify_gp2_elastic(result)
    terms = KC.occupation_reduced_terms(result.occupation)
    @test length(terms) == 9
    @test isempty(KC.reduced_blocked_terms(result.reduced))
    @test isempty(KC.reduced_trotter_terms(result.reduced))

    first_sector = first(keys(terms))
    basis = KC.momentum_basis(first_sector)
    external = KC.external_wigner_momentum(first_sector)
    external_index = only(i for (i, variable) in enumerate(basis) if variable == external)
    loop_indices = Int[i for i in eachindex(basis.variables) if i != external_index]
    @test length(loop_indices) == 2

    k = KC.basis_momentum(basis, external_index)
    q = KC.basis_momentum(basis, loop_indices[1])
    r = KC.basis_momentum(basis, loop_indices[2])
    p = -k + q + r

    ε(momentum) = KC.EnergyForm(KC.DispersionAtom(second_order_pwave_ψ, momentum))
    shell, shell_factor = KC.energy_shell(ε(k) + ε(p) - ε(q) - ε(r))
    @test shell_factor == 1 // 1
    support = KC.FrequencySupport(
        KC.EnergyShell{Fermion}[shell], KC.PrincipalValueSupport{Fermion}[]
    )

    C = KC.ComplexRationals
    n(momentum) =
        KC.OccupationPolynomial(KC.OccupationAtom(second_order_pwave_ψ, momentum), one(C))
    nk, np, nq, nr = n(k), n(p), n(q), n(r)
    one_n = one(typeof(nk))
    fermi_gain_loss =
        (one_n - nk) * (one_n - np) * nq * nr - nk * np * (one_n - nq) * (one_n - nr)

    relative_in = pwave_axis_polynomial(k) - pwave_axis_polynomial(p)
    relative_out = pwave_axis_polynomial(q) - pwave_axis_polynomial(r)
    expected_weight = (1 // 8) * (relative_in * relative_in) * (relative_out * relative_out)

    Monomial = KC.OccupationMonomial{Fermion}
    Polynomial = KC.MomentumPolynomial{KC.ComplexRationals}
    zero_weight = Polynomial()
    generated = Dict{Monomial,Polynomial}()
    for (sector, occupation) in terms
        @test KC.momentum_basis(sector) == basis
        @test KC.external_wigner_momentum(sector) == external
        @test KC.frequency_support(sector) == support
        kinematic = KC.kinematic_factor(sector)
        for (monomial, coefficient) in occupation
            generated[monomial] = get(generated, monomial, zero_weight) + coefficient * kinematic
        end
    end
    filter!(term -> !iszero(last(term)), generated)

    expected = Dict{Monomial,Polynomial}(
        monomial => coefficient * expected_weight for (monomial, coefficient) in fermi_gain_loss
    )
    filter!(term -> !iszero(last(term)), expected)

    @test Set(keys(generated)) == Set(keys(expected))
    for monomial in keys(expected)
        @test generated[monomial] == expected[monomial]
    end
    @test all(
        sector ->
            length(KC.frequency_support(sector).shells) == 1 &&
            isempty(KC.frequency_support(sector).principal_values),
        keys(KC.collision_kernel_terms(result.kernel)),
    )
    return nothing
end

function log_gp2_kernel(result)
    for (sector, occupation) in KC.collision_kernel_terms(result.kernel)
        @info "fermionic p-wave gp² kernel term" basis = KC.momentum_basis(sector) external = KC.external_wigner_momentum(
            sector
        ) kinematic = KC.kinematic_factor(sector) support = KC.frequency_support(sector) occupation
    end
    return nothing
end

@testset "research: spinless-fermion p-wave second-order canonical compiler census" begin
    Lg = fermionic_pwave_elastic_lagrangian()
    Lγ = fermionic_pwave_loss_lagrangian_second_order()
    L = Lg + Lγ

    gp = KC.ParameterMonomial(:gp)
    γp = KC.ParameterMonomial(:γp)
    parameters_expected = (gp^2, gp * γp, γp^2)

    G = DressedPropagator(L, Val(2), Val(5); simplify=true, preserve_regularisation=true)
    @test Set(KC.parameters(G)) == Set(parameters_expected)

    for parameter in parameters_expected
        result = generated_second_order_pipeline(G, parameter)
        @test all(
            object -> KC.statistics(object) === Fermion,
            (
                result.GF,
                result.ΣF,
                result.ΣW,
                result.kinetic,
                result.off_shell,
                result.spectral,
                result.canonical,
                result.reduced,
                result.occupation,
                result.quotient,
                result.kernel,
            ),
        )
        @test KC.target_family(result.spectral) === second_order_pwave_ψ
        @test KC.parameters(result.spectral) == parameter
        @test KC.parameters(result.canonical) == parameter
        @test KC.parameters(result.reduced) == parameter
        @test KC.parameters(result.kernel) == parameter

        census = research_second_order_census(result)
        @info "fermionic p-wave second-order canonical compiler census" parameter census
        if parameter == gp^2
            certify_gp2_elastic(result)
            log_gp2_kernel(result)
        end

        # The research branch may retain genuine blocked/singular support, but every finite
        # strict-QP contribution must traverse the common occupation/loop/kernel path.
        @test census.kernel == census.quotient
        @test census.quotient <= census.occupation
        @test census.occupation <= census.regular
    end
end
