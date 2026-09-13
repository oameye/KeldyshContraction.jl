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
    Pplus_minus, Pplus_plus, Pminus_plus, Pplus_dagger, Pminus_dagger = pwave_branch_pairs(;
        shifted=true
    )

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

function generated_second_order_spectral(G, parameter)
    GF = fourier_transform(G[parameter])
    ΣF = SelfEnergy(GF)
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    kinetic = kinetic_expression(ΣW)
    off_shell = off_shell_collision_expression(kinetic)
    spectral = spectral_dispersive_collision(off_shell)
    return (; GF, ΣF, ΣW, kinetic, off_shell, spectral)
end

function research_term_signature(part, term, coefficient)
    dependency = KC.analyze_spectral_dependencies(term, second_order_pwave_ψ)
    return (;
        part,
        topology=Tuple(KC.topology(term)),
        kinds=Tuple(KC.spectral_dispersive_kinds(term)),
        rank=KC.constraint_rank(dependency),
        count=KC.constraint_count(dependency),
        dependent=KC.has_dependent_shell_support(dependency),
        loops=KC.loop_frequency_count(term),
        shifts=Tuple(regularisation_shift(line) for line in kinetic_lines(term.carrier)),
        coefficient,
        kinematic=KC.kinematic_factor(term),
    )
end

function probe_frequency_reduction(collision)
    successes = NamedTuple[]
    failures = NamedTuple[]
    for (part, expression) in (
        (:offset, collision_offset(collision)),
        (:distribution, collision_distribution_coefficient(collision)),
    )
        for (term, coefficient) in expression
            signature = research_term_signature(part, term, coefficient)
            try
                result = KC.reduce_frequency_term(term, second_order_pwave_ψ)
                push!(
                    successes,
                    merge(
                        signature,
                        (;
                            regular=length(result.regular),
                            dependent_result=length(result.dependent),
                            causal=length(result.causal),
                            trotter=length(result.trotter),
                            zero=isempty(result),
                        ),
                    ),
                )
            catch error
                push!(
                    failures,
                    merge(
                        signature,
                        (;
                            error_type=typeof(error), error_message=sprint(showerror, error)
                        ),
                    ),
                )
            end
        end
    end
    return successes, failures
end

function reduction_census(records)
    census = Dict{Any,Int}()
    for record in records
        key = (
            record.part,
            record.topology,
            record.kinds,
            record.rank,
            record.count,
            record.dependent,
            record.loops,
            record.regular,
            record.dependent_result,
            record.causal,
            record.trotter,
            record.zero,
        )
        census[key] = get(census, key, 0) + 1
    end
    return census
end

@testset "research: complete spinless-fermion p-wave second-order census" begin
    Lg = fermionic_pwave_elastic_lagrangian()
    Lγ = fermionic_pwave_loss_lagrangian_second_order()
    L = Lg + Lγ

    gp = KC.ParameterMonomial(:gp)
    γp = KC.ParameterMonomial(:γp)
    gp² = gp^2
    gpγp = gp * γp
    γp² = γp^2

    G = DressedPropagator(L, Val(2), Val(5); simplify=true, preserve_regularisation=true)
    @test Set(KC.parameters(G)) == Set((gp², gpγp, γp²))

    for parameter in (gp², gpγp, γp²)
        result = generated_second_order_spectral(G, parameter)
        @test all(
            object -> KC.statistics(object) === Fermion,
            (
                result.GF,
                result.ΣF,
                result.ΣW,
                result.kinetic,
                result.off_shell,
                result.spectral,
            ),
        )
        @test target_family(result.spectral) === second_order_pwave_ψ
        @test parameters(result.spectral) == parameter

        successes, failures = probe_frequency_reduction(result.spectral)
        @info "fermionic p-wave second-order frequency census" parameter spectral_terms = length(
            successes
        ) + length(failures) reduction_census = reduction_census(successes) failures
    end
end
