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

function certify_compact_second_order_display(result)
    terms = KC.collision_kernel_terms(result.kernel)
    expanded = KC._physical_collision_string(terms, true)
    compact = KC._compact_physical_collision_string(terms, true)

    @test ncodeunits(compact) < ncodeunits(expanded)
    @test length(collect(eachmatch(r"_\{x\}\^2", compact))) >= 2
    @test !occursin(r"q_\{\d+\}_\{", compact)
    return compact
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

    # The physical Boltzmann form is an equivalence class under dummy-loop relabeling.
    # Build it independently, split only the kinematic polynomial into exact monomials,
    # then send it through the same safe quotient before comparing final kernels.
    Sector = typeof(first_sector)
    Occupation = typeof(fermi_gain_loss)
    expected_terms = Dict{Sector,Occupation}()
    zero_occupation = zero(fermi_gain_loss)
    for (monomial, coefficient) in expected_weight
        kinematic = KC.MomentumPolynomial(monomial, one(C))
        sector = KC.ReducedCollisionSector(
            KC.parameters(result.occupation), basis, external, kinematic, support
        )
        contribution = coefficient * fermi_gain_loss
        expected_terms[sector] = get(expected_terms, sector, zero_occupation) + contribution
    end

    expected_occupation = typeof(result.occupation)(
        expected_terms,
        KC.target_family(result.occupation),
        KC.parameters(result.occupation),
        KC.wigner_context(result.occupation),
    )
    expected_kernel = KC.collision_kernel(KC.quotient_loop_momenta(expected_occupation))

    @test KC.collision_kernel_terms(result.kernel) ==
        KC.collision_kernel_terms(expected_kernel)
    @test all(
        sector ->
            length(KC.frequency_support(sector).shells) == 1 &&
            isempty(KC.frequency_support(sector).principal_values),
        keys(KC.collision_kernel_terms(result.kernel)),
    )

    compact = certify_compact_second_order_display(result)
    @test occursin("\\left(1-n_{", compact)
    return nothing
end

function certify_gamma2_loss(result)
    terms = KC.occupation_reduced_terms(result.occupation)
    @test length(terms) == 9
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

    # The trace-preserving dissipative pair cut is fixed directly by the fermionic
    # greater/lesser pair weights, Π>=(1-n₁)(1-n₂) and Π<=n₁n₂.
    pair_cut =
        (one_n - nk) * (one_n - np) * nq * nr + nk * np * (one_n - nq) * (one_n - nr) -
        2 * nk * np * nq * nr

    relative_in = pwave_axis_polynomial(k) - pwave_axis_polynomial(p)
    relative_out = pwave_axis_polynomial(q) - pwave_axis_polynomial(r)
    expected_weight = (1 // 8) * (relative_in * relative_in) * (relative_out * relative_out)

    Sector = typeof(first_sector)
    Occupation = typeof(pair_cut)
    expected_terms = Dict{Sector,Occupation}()
    zero_occupation = zero(pair_cut)
    for (monomial, coefficient) in expected_weight
        kinematic = KC.MomentumPolynomial(monomial, one(C))
        sector = KC.ReducedCollisionSector(
            KC.parameters(result.occupation), basis, external, kinematic, support
        )
        contribution = coefficient * pair_cut
        expected_terms[sector] = get(expected_terms, sector, zero_occupation) + contribution
    end

    expected_occupation = typeof(result.occupation)(
        expected_terms,
        KC.target_family(result.occupation),
        KC.parameters(result.occupation),
        KC.wigner_context(result.occupation),
    )
    expected_kernel = KC.collision_kernel(KC.quotient_loop_momenta(expected_occupation))

    @test KC.collision_kernel_terms(result.kernel) ==
        KC.collision_kernel_terms(expected_kernel)
    @test all(
        sector ->
            length(KC.frequency_support(sector).shells) == 1 &&
            isempty(KC.frequency_support(sector).principal_values),
        keys(KC.collision_kernel_terms(result.kernel)),
    )

    compact = certify_compact_second_order_display(result)
    @test occursin("\\left(1-n_{", compact)

    # The singular branch remains upstream. Every generated blocker is a genuine pinch,
    # and all 84 contributions carry the same exact canonical affine singular geometry.
    blocked = collect(values(KC.reduced_blocked_terms(result.reduced)))
    @test length(blocked) == 84
    @test all(
        contribution ->
            KC.trotter_frequency_blocker_kind(contribution) === KC.CausalFrequencyPinch,
        blocked,
    )
    blocked_supports = [
        KC.trotter_frequency_blocked_supports(contribution) for contribution in blocked
    ]
    reference_supports = first(blocked_supports)
    @test !isempty(reference_supports)
    @test all(KC.has_affine_singular_support, reference_supports)
    @test all(supports -> supports == reference_supports, blocked_supports)
    @info "fermionic p-wave gamma2 blocked causal census" blocker_count = length(blocked) support_ranks = [
        KC.affine_support_rank(support) for support in reference_supports
    ] support_counts = [KC.affine_support_count(support) for support in reference_supports]
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
        parameter == gp^2 && certify_gp2_elastic(result)
        parameter == γp^2 && certify_gamma2_loss(result)

        # The research branch may retain genuine blocked/singular support, but every finite
        # strict-QP contribution must traverse the common occupation/loop/kernel path.
        @test census.kernel == census.quotient
        @test census.quotient <= census.occupation
        @test census.occupation <= census.regular
    end
end