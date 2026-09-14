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
        key = repr([
            affine_support_summary_gamma2(support) for support in supports
        ])
        support_counts[key] = get(support_counts, key, 0) + 1
    end
    @info "PWAVE_GAMMA2_SUPPORT_COUNTS" support_counts

    @test length(regular) == 9
    @test length(blocked) == 84
end
