using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields diagnostic_gamma2_ψ::Fermion

function diagnostic_gamma2_loss_lagrangian()
    ψ₁ = diagnostic_gamma2_ψ[One]
    ψ₂ = diagnostic_gamma2_ψ[Two]
    bψ₁ = bar(ψ₁)
    bψ₂ = bar(ψ₂)
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus

    ψplus_minus = ψ₁(minus) + ψ₂(minus)
    ∂ψplus_minus = partial(ψ₁(minus), :x) + partial(ψ₂(minus), :x)
    ψplus_plus = ψ₁(plus) + ψ₂(plus)
    ∂ψplus_plus = partial(ψ₁(plus), :x) + partial(ψ₂(plus), :x)
    ψminus_plus = ψ₁(plus) - ψ₂(plus)
    ∂ψminus_plus = partial(ψ₁(plus), :x) - partial(ψ₂(plus), :x)

    bψplus = bψ₁ + bψ₂
    bψminus = bψ₂ - bψ₁
    ∂bψplus = partial(bψ₁, :x) + partial(bψ₂, :x)
    ∂bψminus = partial(bψ₂, :x) - partial(bψ₁, :x)

    Pplus_minus = ψplus_minus * ∂ψplus_minus
    Pplus_plus = ψplus_plus * ∂ψplus_plus
    Pminus_plus = ψminus_plus * ∂ψminus_plus
    Pplus_dagger = ∂bψplus * bψplus
    Pminus_dagger = ∂bψminus * bψminus

    loss =
        (1 // 8) *
        im *
        (
            (Pplus_dagger - Pminus_dagger) * Pplus_minus -
            (Pplus_plus - Pminus_plus) * Pminus_dagger
        )
    return InteractionLagrangian(loss, :γp)
end

momentum_tuple_gamma2(momentum) = Tuple(momentum.coefficients)

function gamma2_axis(momentum)
    C = KC.ComplexRationals
    component = KC.MomentumComponent(momentum, :x)
    monomial = KC.MomentumMonomial(KC.MomentumComponent[component])
    return KC.MomentumPolynomial(monomial, one(C))
end

function kinematic_summary_gamma2(polynomial)
    return [
        (
            coefficient=coefficient,
            factors=Tuple(
                (
                    axis=factor.axis,
                    momentum=momentum_tuple_gamma2(factor.momentum),
                ) for factor in monomial
            ),
        ) for (monomial, coefficient) in polynomial
    ]
end

function occupation_summary_gamma2(polynomial)
    return [
        (
            coefficient=coefficient,
            factors=Tuple(momentum_tuple_gamma2(atom.momentum) for atom in monomial),
        ) for (monomial, coefficient) in polynomial
    ]
end

function energy_summary_gamma2(form)
    return [
        (coefficient=coefficient, momentum=momentum_tuple_gamma2(atom.momentum)) for
        (atom, coefficient) in KC.energy_terms(form)
    ]
end

function support_summary_gamma2(support)
    return (
        shells=[energy_summary_gamma2(shell.energy) for shell in support.shells],
        principal_values=[energy_summary_gamma2(pv.energy) for pv in support.principal_values],
    )
end

function affine_support_summary_gamma2(support)
    return (
        independent=[
            (
                loop=Tuple(constraint.loop_coefficients),
                energy=energy_summary_gamma2(constraint.energy),
            ) for constraint in support.independent_constraints
        ],
        dependencies=[Tuple(row) for row in support.dependency_rows],
    )
end

function expected_gamma2_kernel(occupation)
    terms = KC.occupation_reduced_terms(occupation)
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

    ε(momentum) = KC.EnergyForm(KC.DispersionAtom(diagnostic_gamma2_ψ, momentum))
    shell, shell_factor = KC.energy_shell(ε(k) + ε(p) - ε(q) - ε(r))
    @test shell_factor == 1 // 1
    support = KC.FrequencySupport(
        KC.EnergyShell{Fermion}[shell], KC.PrincipalValueSupport{Fermion}[]
    )

    C = KC.ComplexRationals
    n(momentum) =
        KC.OccupationPolynomial(KC.OccupationAtom(diagnostic_gamma2_ψ, momentum), one(C))
    nk, np, nq, nr = n(k), n(p), n(q), n(r)
    one_n = one(typeof(nk))

    pair_cut =
        (one_n - nk) * (one_n - np) * nq * nr +
        nk * np * (one_n - nq) * (one_n - nr) -
        2 * nk * np * nq * nr

    relative_in = gamma2_axis(k) - gamma2_axis(p)
    relative_out = gamma2_axis(q) - gamma2_axis(r)
    expected_weight =
        (1 // 8) * (relative_in * relative_in) * (relative_out * relative_out)

    Sector = typeof(first_sector)
    Occupation = typeof(pair_cut)
    expected_terms = Dict{Sector,Occupation}()
    zero_occupation = zero(pair_cut)
    for (monomial, coefficient) in expected_weight
        kinematic = KC.MomentumPolynomial(monomial, one(C))
        sector = KC.ReducedCollisionSector(
            KC.ParameterMonomial(:γp)^2, basis, external, kinematic, support
        )
        contribution = coefficient * pair_cut
        expected_terms[sector] = get(expected_terms, sector, zero_occupation) + contribution
    end

    expected_occupation = typeof(occupation)(
        expected_terms,
        KC.target_family(occupation),
        KC.parameters(occupation),
        KC.wigner_context(occupation),
    )
    return KC.collision_kernel(KC.quotient_loop_momenta(expected_occupation))
end

@testset "diagnostic: spinless-fermion p-wave gamma2 sector" begin
    L = diagnostic_gamma2_loss_lagrangian()
    G = DressedPropagator(L, Val(2), Val(5); simplify=true, preserve_regularisation=true)
    @test KC.parameters(G) == KC.ParameterMonomial(:γp)^2
    GF = fourier_transform(G)
    ΣW = wigner_transform(SelfEnergy(GF); gradient_order=Val(0))
    spectral = spectral_dispersive_collision(
        off_shell_collision_expression(kinetic_expression(ΣW))
    )
    canonical = KC.canonical_frequency_collision(spectral)
    reduced = KC.reduce_frequency_collision(canonical)
    occupation = KC.occupation_reduced_expression(reduced)
    quotient = KC.quotient_loop_momenta(occupation)
    kernel = KC.collision_kernel(quotient)

    regular = [
        (
            basis=length(KC.momentum_basis(sector)),
            external=KC.external_wigner_momentum(sector).index,
            kinematic=kinematic_summary_gamma2(KC.kinematic_factor(sector)),
            support=support_summary_gamma2(KC.frequency_support(sector)),
            occupation=occupation_summary_gamma2(polynomial),
        ) for (sector, polynomial) in KC.loop_quotient_terms(quotient)
    ]
    sort!(regular; by=repr)
    @info "PWAVE_GAMMA2_REGULAR" regular

    expected = expected_gamma2_kernel(occupation)
    @test KC.collision_kernel_terms(kernel) == KC.collision_kernel_terms(expected)

    blocked = collect(values(KC.reduced_blocked_terms(reduced)))
    kind_counts = Dict{String,Int}()
    for contribution in blocked
        key = string(KC.trotter_frequency_blocker_kind(contribution))
        kind_counts[key] = get(kind_counts, key, 0) + 1
    end
    @info "PWAVE_GAMMA2_BLOCKER_COUNTS" kind_counts

    support_counts = Dict{String,Int}()
    for contribution in blocked
        supports = KC.trotter_frequency_blocked_supports(contribution)
        key = repr([affine_support_summary_gamma2(support) for support in supports])
        support_counts[key] = get(support_counts, key, 0) + 1
    end
    @info "PWAVE_GAMMA2_SUPPORT_COUNTS" support_counts

    @test length(regular) == 9
    @test length(blocked) == 84
end
