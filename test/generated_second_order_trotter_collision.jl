using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields canonical_trotter_loss2_ϕ::Boson

function generated_second_order_canonical_trotter_collision()
    c = canonical_trotter_loss2_ϕ[Classical]
    q = canonical_trotter_loss2_ϕ[Quantum]
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

    L = @inferred InteractionLagrangian(loss, :γ)
    G = @inferred DressedPropagator(
        L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false
    )
    GF = @inferred fourier_transform(G)
    ΣF = @inferred SelfEnergy(GF)
    ΣW = @inferred wigner_transform(ΣF; gradient_order=Val(0))
    KΣ = @inferred kinetic_expression(ΣW)
    off_shell = @inferred off_shell_collision_expression(KΣ)
    spectral = @inferred spectral_dispersive_collision(off_shell)
    return @inferred KC.canonical_frequency_collision(spectral)
end

function second_order_trotter_support_signature(supports)
    return Tuple(
        (
            KC.affine_support_rank(support),
            KC.affine_support_count(support),
            Tuple(Tuple(row) for row in support.dependency_rows),
        ) for support in supports
    )
end

function second_order_boundary_has_no_active_frequency(expression)
    return all(KC.causal_frequency_terms(expression)) do term
        all(KC.causal_frequency_denominators(term)) do denominator
            return all(iszero, denominator.loop_coefficients)
        end
    end
end

@testset "generated γ² grouped canonical Trotter reduction" begin
    collision = generated_second_order_canonical_trotter_collision()
    @test parameters(collision) == KC.ParameterMonomial(:γ)^2
    @test target_family(collision) === canonical_trotter_loss2_ϕ

    shifted = KC.shifted_frequency_terms(collision)
    @test length(shifted) == 32

    isolated = 0
    for term in shifted
        state = @inferred KC.TrotterFrequencyState(term)
        witness = @inferred KC.trotter_frequency_witness(state)
        if KC.has_isolated_trotter_frequency(witness)
            isolated += 1
            line = KC.trotter_frequency_line(state, witness)
            @test statistical_weight(line) === NoStatisticalWeight
        end
        reduced = @inferred KC.eliminate_isolated_trotter_frequencies(state)
        @test !KC.has_unresolved_trotter_frequency(reduced)
    end
    @test isolated == 32

    grouped = @inferred KC.reduce_shifted_trotter_frequencies(collision)
    @test isempty(KC.unresolved_trotter_states(grouped))

    expressions = KC.grouped_trotter_expressions(grouped)
    constants = KC.grouped_trotter_constants(grouped)
    @test length(expressions) == 4
    @test isempty(constants)
    @test all(!isempty, values(expressions))

    reversed_collision = typeof(collision)(
        copy(KC.canonical_frequency_expressions(collision)),
        reverse(shifted),
        target_family(collision),
        parameters(collision),
        KC.wigner_context(collision),
    )
    reversed_grouped = @inferred KC.reduce_shifted_trotter_frequencies(reversed_collision)
    @test isempty(KC.unresolved_trotter_states(reversed_grouped))
    @test KC.grouped_trotter_expressions(reversed_grouped) == expressions
    @test KC.grouped_trotter_constants(reversed_grouped) == constants

    # Equal-depth lower stages are selected canonically rather than through Dict iteration order.
    sector, expression = first(KC.canonical_frequency_expressions(collision))
    canonical = KC._complex_causal_frequency_expression(expression)
    source_key = KC.TrotterFrequencyGroupKey(sector, [1, 2])
    key_1 = KC.TrotterFrequencyGroupKey(sector, [1])
    key_2 = KC.TrotterFrequencyGroupKey(sector, [2])
    stages_21 = Dict(key_2 => canonical, key_1 => canonical)
    stages_12 = Dict(key_1 => canonical, key_2 => canonical)
    found_21, subset_21 = @inferred KC._nearest_existing_trotter_subset(
        stages_21, source_key
    )
    found_12, subset_12 = @inferred KC._nearest_existing_trotter_subset(
        stages_12, source_key
    )
    @test found_21
    @test found_12
    @test KC.active_trotter_frequencies(subset_21) == [1]
    @test KC.active_trotter_frequencies(subset_12) == [1]
    @test KC._next_trotter_frequency(stages_21, source_key) == 2
    @test KC._next_trotter_frequency(stages_12, source_key) == 2

    reduction = @inferred KC.reduce_canonical_trotter_frequencies(collision)
    @test isempty(KC.unresolved_trotter_states(reduction))
    blockers = KC.blocked_trotter_contributions(reduction)
    reduced_constants = KC.canonical_trotter_constants(reduction)
    boundary = KC.canonical_trotter_boundary_expressions(reduction)

    @test isempty(reduced_constants)
    # The eight regular sunset reductions are all retained. Four were previously quarantined
    # inside pinched complete sectors, which was the γ² regression caught by the raw oracle.
    @test length(boundary) == 8
    @test all(second_order_boundary_has_no_active_frequency, values(boundary))
    @test length(blockers) == 8

    active_signatures = Tuple(
        Tuple(KC.active_trotter_frequencies(blocked)) for blocked in values(blockers)
    )
    @test count(==((1,)), active_signatures) == 4
    @test count(==((1, 2)), active_signatures) == 4

    expected_support = ((1, 2, ((1 // 1,),)),)
    for blocked in values(blockers)
        @test KC.trotter_frequency_blocker_kind(blocked) === KC.CausalFrequencyPinch
        @test KC.trotter_frequency_blocked_index(blocked) == 1
        @test second_order_trotter_support_signature(
            KC.trotter_frequency_blocked_supports(blocked)
        ) == expected_support
    end

    reversed_reduction = @inferred KC.reduce_canonical_trotter_frequencies(
        reversed_collision
    )
    @test isempty(KC.unresolved_trotter_states(reversed_reduction))
    @test KC.canonical_trotter_constants(reversed_reduction) == reduced_constants
    @test KC.canonical_trotter_boundary_expressions(reversed_reduction) == boundary

    reversed_blockers = KC.blocked_trotter_contributions(reversed_reduction)
    @test Set(keys(reversed_blockers)) == Set(keys(blockers))
    for (key, blocked) in blockers
        reversed_blocked = reversed_blockers[key]
        @test KC.trotter_frequency_blocker_kind(reversed_blocked) ===
            KC.trotter_frequency_blocker_kind(blocked)
        @test KC.trotter_frequency_blocked_index(reversed_blocked) ==
            KC.trotter_frequency_blocked_index(blocked)
        @test KC.trotter_frequency_blocked_expression(reversed_blocked) ==
            KC.trotter_frequency_blocked_expression(blocked)
        @test second_order_trotter_support_signature(
            KC.trotter_frequency_blocked_supports(reversed_blocked)
        ) == second_order_trotter_support_signature(
            KC.trotter_frequency_blocked_supports(blocked)
        )
    end
end
