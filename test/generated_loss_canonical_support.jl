using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields canonical_loss_support_ϕ::Boson

function generated_loss_canonical_frequency_collision()
    c = canonical_loss_support_ϕ[Classical]
    q = canonical_loss_support_ϕ[Quantum]
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

function affine_support_signature(support)
    independent = Tuple(
        (Tuple(constraint.loop_coefficients), Tuple(KC.energy_terms(constraint.energy))) for
        constraint in support.independent_constraints
    )
    dependencies = Tuple(Tuple(row) for row in support.dependency_rows)
    return (;
        rank=KC.affine_support_rank(support),
        count=KC.affine_support_count(support),
        independent,
        dependencies,
    )
end

@testset "generated γ² canonical singular support" begin
    collision = generated_loss_canonical_frequency_collision()
    @test parameters(collision) == KC.ParameterMonomial(:γ)^2
    @test target_family(collision) === canonical_loss_support_ϕ

    shifted = KC.shifted_frequency_terms(collision)
    expressions = KC.canonical_frequency_expressions(collision)
    @test !isempty(shifted)
    @test !isempty(expressions)

    forward_plan = KC.CausalFrequencyReductionPlan((1, 2))
    reverse_plan = KC.CausalFrequencyReductionPlan((2, 1))
    census = NamedTuple[]

    for (sector, expression) in expressions
        nfrequencies = length(KC.momentum_basis(sector)) - 1
        @test nfrequencies == 2

        forward = KC.reduce_causal_frequency_with_support(expression, forward_plan)
        reverse = KC.reduce_causal_frequency_with_support(expression, reverse_plan)
        forward_kind = KC.causal_frequency_support_kind(forward)
        reverse_kind = KC.causal_frequency_support_kind(reverse)
        forward_supports = KC.causal_frequency_singular_supports(forward)
        reverse_supports = KC.causal_frequency_singular_supports(reverse)
        forward_signatures = Tuple(affine_support_signature(s) for s in forward_supports)
        reverse_signatures = Tuple(affine_support_signature(s) for s in reverse_supports)

        @test forward_kind === reverse_kind
        @test forward_kind === KC.CausalFrequencyIntegrated ||
            forward_kind === KC.CausalFrequencyPinch

        if forward_kind === KC.CausalFrequencyIntegrated
            @test KC.causal_frequency_support_expression(forward) ==
                KC.causal_frequency_support_expression(reverse)
            @test isempty(forward_supports)
            @test isempty(reverse_supports)
        else
            @test !isempty(forward_supports)
            @test forward_signatures == reverse_signatures
        end

        statistical = Tuple(
            Tuple(KC.momentum(atom).coefficients) for
            atom in KC.statistical_monomial(sector)
        )
        push!(
            census,
            (;
                statistical,
                source_terms=length(KC.causal_frequency_terms(expression)),
                forward_kind,
                reverse_kind,
                forward_blocked=KC.causal_frequency_support_blocked_frequency(forward),
                reverse_blocked=KC.causal_frequency_support_blocked_frequency(reverse),
                forward_supports=forward_signatures,
                reverse_supports=reverse_signatures,
            ),
        )
    end

    pinch = filter(record -> record.forward_kind === KC.CausalFrequencyPinch, census)
    regular = filter(record -> record.forward_kind === KC.CausalFrequencyIntegrated, census)

    # Complete statistical-sector assembly and every regular contour integration reduce the
    # unshifted γ² sector to seven regular sectors and four genuine pinch sectors. In each pinch,
    # the independent third shell is a coloop and therefore remains ordinary residual support;
    # the unresolved singular subsystem is the rank-one repeated shell A(q)^2.
    @test length(pinch) == 4
    @test length(regular) == 7
    @test all(pinch) do record
        length(record.forward_supports) == 1 || return false
        support = only(record.forward_supports)
        return support.rank == 1 &&
               support.count == 2 &&
               support.dependencies == ((1 // 1,),) &&
               record.forward_supports == record.reverse_supports
    end

    @info "generated γ² canonical support" shifted_count=length(shifted) regular_count=length(
        regular
    ) pinch_count=length(pinch)
end
