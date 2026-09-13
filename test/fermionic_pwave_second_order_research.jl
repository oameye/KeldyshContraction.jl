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
    ∂ψplus = shifted ?
        partial(ψ₁(minus), :x) + partial(ψ₂(minus), :x) :
        partial(ψ₁, :x) + partial(ψ₂, :x)

    ψplus_forward = shifted ? ψ₁(plus) + ψ₂(plus) : ψ₁ + ψ₂
    ∂ψplus_forward = shifted ?
        partial(ψ₁(plus), :x) + partial(ψ₂(plus), :x) :
        partial(ψ₁, :x) + partial(ψ₂, :x)

    ψminus = shifted ? ψ₁(plus) - ψ₂(plus) : ψ₁ - ψ₂
    ∂ψminus = shifted ?
        partial(ψ₁(plus), :x) - partial(ψ₂(plus), :x) :
        partial(ψ₁, :x) - partial(ψ₂, :x)

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
    Pplus, _, Pminus, Pplus_dagger, Pminus_dagger = pwave_branch_pairs(
        shifted=false
    )

    # Physical convention:
    #
    #   H_int = (g_p/2) P† P,       P = ψ ∂x ψ,
    #   S_g   = -(g_p/2) ∫ (P+†P+ - P-†P-).
    #
    # Each branch pair contains two LO factors 1/√2. The unnormalised LO sums above
    # therefore leave an exact coefficient -1/8 after the branch rotation.
    elastic = -(1 // 8) * (Pplus_dagger * Pplus - Pminus_dagger * Pminus)
    return InteractionLagrangian(elastic, :gp)
end

function fermionic_pwave_loss_lagrangian_second_order()
    Pplus_minus, Pplus_plus, Pminus_plus, Pplus_dagger, Pminus_dagger =
        pwave_branch_pairs(; shifted=true)

    # Same finite-Trotter Lindblad convention as the frozen O(γp) oracle in #312.
    loss =
        (1 // 8) *
        im *
        (
            (Pplus_dagger - Pminus_dagger) * Pplus_minus -
            (Pplus_plus - Pminus_plus) * Pminus_dagger
        )
    return InteractionLagrangian(loss, :γp)
end

function component_topology_census(Σ)
    census(component) = Dict(
        Tuple(topology) => length(diagrams) for (topology, diagrams) in topologies(component)
    )
    return (
        keldysh=census(KC.keldysh_component(Σ)),
        retarded=census(KC.retarded_component(Σ)),
        advanced=census(KC.advanced_component(Σ)),
    )
end

function generated_second_order_sector(G, parameter)
    GF = fourier_transform(G[parameter])
    ΣF = SelfEnergy(GF)
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    kinetic = kinetic_expression(ΣW)
    off_shell = off_shell_collision_expression(kinetic)
    spectral = spectral_dispersive_collision(off_shell)
    reduced = reduce_frequency_collision(spectral)
    occupation = occupation_reduced_expression(reduced)
    quotient = quotient_loop_momenta(occupation)
    kernel = collision_kernel(quotient)
    return (; GF, ΣF, ΣW, kinetic, off_shell, spectral, reduced, occupation, quotient, kernel)
end

@testset "research: complete spinless-fermion p-wave second-order census" begin
    Lg = fermionic_pwave_elastic_lagrangian()
    Lγ = fermionic_pwave_loss_lagrangian_second_order()
    L = Lg + Lγ

    @test KC.statistics(L) === Fermion
    @test KC.target_family(L) === second_order_pwave_ψ

    gp = KC.ParameterMonomial(:gp)
    γp = KC.ParameterMonomial(:γp)
    gp² = gp^2
    gpγp = gp * γp
    γp² = γp^2

    G = DressedPropagator(
        L, Val(2), Val(5); simplify=true, preserve_regularisation=true
    )
    @test Set(KC.parameters(G)) == Set((gp², gpγp, γp²))

    for parameter in (gp², gpγp, γp²)
        result = generated_second_order_sector(G, parameter)
        @test all(
            object -> KC.statistics(object) === Fermion,
            (
                result.GF,
                result.ΣF,
                result.ΣW,
                result.kinetic,
                result.off_shell,
                result.spectral,
                result.reduced,
                result.occupation,
                result.quotient,
                result.kernel,
            ),
        )
        @test KC.parameters(result.kernel) == parameter

        regular = KC.reduced_regular_terms(result.reduced)
        dependent = KC.reduced_dependent_terms(result.reduced)
        causal = KC.reduced_causal_terms(result.reduced)
        trotter = KC.reduced_trotter_terms(result.reduced)

        @info "fermionic p-wave second-order sector" parameter topology_census = component_topology_census(
            result.ΣF
        ) regular_count = length(regular) dependent_count = length(dependent) causal_count = length(
            causal
        ) trotter_count = length(trotter) occupation_terms = collect(
            KC.occupation_reduced_terms(result.occupation)
        ) kernel_terms = collect(KC.collision_kernel_terms(result.kernel))
    end
end
