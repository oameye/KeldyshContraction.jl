using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields mixed_pwave_ψ::Fermion

function mixed_pwave_branch_pairs(; shifted::Bool)
    ψ₁ = mixed_pwave_ψ[One]
    ψ₂ = mixed_pwave_ψ[Two]
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

function mixed_pwave_lagrangian()
    Pplus, _, Pminus, Pplus_dagger, Pminus_dagger = mixed_pwave_branch_pairs(;
        shifted=false
    )
    elastic = -(1 // 8) * (Pplus_dagger * Pplus - Pminus_dagger * Pminus)

    Pplus_minus, Pplus_plus, Pminus_plus, Pplus_dagger_shifted, Pminus_dagger_shifted = mixed_pwave_branch_pairs(;
        shifted=true
    )
    loss =
        (1 // 8) *
        im *
        (
            (Pplus_dagger_shifted - Pminus_dagger_shifted) * Pplus_minus -
            (Pplus_plus - Pminus_plus) * Pminus_dagger_shifted
        )

    return InteractionLagrangian(elastic, :gp) + InteractionLagrangian(loss, :γp)
end

function mixed_pwave_pipeline()
    L = mixed_pwave_lagrangian()
    G = DressedPropagator(L, Val(2), Val(5); simplify=true, preserve_regularisation=true)
    parameter = KC.ParameterMonomial(:gp) * KC.ParameterMonomial(:γp)
    GF = fourier_transform(G[parameter])
    ΣW = wigner_transform(SelfEnergy(GF); gradient_order=Val(0))
    spectral = spectral_dispersive_collision(
        off_shell_collision_expression(kinetic_expression(ΣW))
    )
    canonical = KC.canonical_frequency_collision(spectral)
    reduced = KC.reduce_frequency_collision(canonical)
    occupation = KC.occupation_reduced_expression(reduced)
    quotient = KC.quotient_loop_momenta(occupation)
    kernel = KC.collision_kernel(quotient)
    return (; parameter, spectral, canonical, reduced, occupation, quotient, kernel)
end

function mixed_pwave_axis(momentum)
    C = KC.ComplexRationals
    component = KC.MomentumComponent(momentum, :x)
    monomial = KC.MomentumMonomial(KC.MomentumComponent[component])
    return KC.MomentumPolynomial(monomial, one(C))
end

function expected_mixed_pwave_kernel(result)
    terms = KC.occupation_reduced_terms(result.occupation)
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

    ε(momentum) = KC.EnergyForm(KC.DispersionAtom(mixed_pwave_ψ, momentum))
    ΔE = ε(k) + ε(p) - ε(q) - ε(r)
    pv, pv_factor = KC.principal_value_support(ΔE)
    support = KC.FrequencySupport(
        KC.EnergyShell{Fermion}[], KC.PrincipalValueSupport{Fermion}[pv]
    )

    C = KC.ComplexRationals
    n(momentum) =
        KC.OccupationPolynomial(KC.OccupationAtom(mixed_pwave_ψ, momentum), one(C))
    nk, np, nq, nr = n(k), n(p), n(q), n(r)
    one_n = one(typeof(nk))

    # Independent medium pair factor. For a fermionic virtual pair the retarded pair bubble
    # carries 1 - n_q - n_r = (F_q + F_r)/2, not a final-state Boltzmann blocking factor.
    pair_medium = nk * np * (one_n - nq - nr)

    relative_in = mixed_pwave_axis(k) - mixed_pwave_axis(p)
    relative_virtual = mixed_pwave_axis(q) - mixed_pwave_axis(r)

    # H_eff contains λ_p P†P/2 with λ_p = g_p - iγ_p. The unordered intermediate pair
    # contributes λ_p²/8. Since 2 Im(λ_p²/8) = -g_p γ_p/2, the regular mixed term is
    # purely dispersive and has coefficient -1/2 multiplying PV(1/ΔE).
    expected_weight =
        -(1 // 2) * (relative_in * relative_in) * (relative_virtual * relative_virtual)

    Sector = typeof(first_sector)
    Occupation = typeof(pair_medium)
    expected_terms = Dict{Sector,Occupation}()
    zero_occupation = zero(pair_medium)
    for (monomial, coefficient) in expected_weight
        kinematic = KC.MomentumPolynomial(monomial, one(C))
        sector = KC.ReducedCollisionSector(
            result.parameter, basis, external, kinematic, support
        )
        contribution = coefficient * pv_factor * pair_medium
        expected_terms[sector] = get(expected_terms, sector, zero_occupation) + contribution
    end

    expected_occupation = typeof(result.occupation)(
        expected_terms,
        KC.target_family(result.occupation),
        result.parameter,
        KC.wigner_context(result.occupation),
    )
    return KC.collision_kernel(KC.quotient_loop_momenta(expected_occupation))
end

@testset "research: spinless-fermion p-wave mixed gpγp regular kernel" begin
    result = mixed_pwave_pipeline()

    @test length(KC.occupation_reduced_terms(result.occupation)) == 9
    @test length(KC.reduced_blocked_terms(result.reduced)) == 30
    @test isempty(KC.reduced_trotter_terms(result.reduced))

    # At O(g_p γ_p) the on-shell Born term vanishes because
    # |g_p - iγ_p|² = g_p² + γ_p². Any finite mixed contribution is therefore PV-only.
    @test all(
        sector ->
            isempty(KC.frequency_support(sector).shells) &&
            length(KC.frequency_support(sector).principal_values) == 1,
        keys(KC.occupation_reduced_terms(result.occupation)),
    )

    expected = expected_mixed_pwave_kernel(result)
    @test KC.collision_kernel_terms(result.kernel) == KC.collision_kernel_terms(expected)

    blocker_counts = Dict{KC.CausalFrequencyIntegrationKind,Int}()
    for blocked in values(KC.reduced_blocked_terms(result.reduced))
        kind = KC.trotter_frequency_blocker_kind(blocked)
        blocker_counts[kind] = get(blocker_counts, kind, 0) + 1
    end
    @info "fermionic p-wave gpγp blocked causal census" blocker_counts
end
