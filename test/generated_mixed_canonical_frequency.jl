using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields canonical_mixed_ϕ::Boson

function generated_mixed_canonical_collision()
    c = canonical_mixed_ϕ[Classical]
    q = canonical_mixed_ϕ[Quantum]
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus

    elastic = -(1 // 2) * ((c^2 + q^2) * bar(c) * bar(q) + c * q * (bar(c)^2 + bar(q)^2))
    loss =
        (1 // 2) *
        im *
        (
            bar(c) * bar(q) * (c(minus)^2 + q(minus)^2) -
            c(plus) * q(plus) * (bar(c)^2 + bar(q)^2) +
            2 * bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
        )

    L = @inferred InteractionLagrangian(elastic, :g) + InteractionLagrangian(loss, :γ)
    G = @inferred DressedPropagator(
        L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false
    )
    mixed_parameter = KC.ParameterMonomial(:g) * KC.ParameterMonomial(:γ)
    GF = @inferred fourier_transform(G[mixed_parameter])
    ΣF = @inferred SelfEnergy(GF)
    ΣW = @inferred wigner_transform(ΣF; gradient_order=Val(0))
    KΣ = @inferred kinetic_expression(ΣW)
    off_shell = @inferred off_shell_collision_expression(KΣ)
    collision = @inferred spectral_dispersive_collision(off_shell)
    return @inferred KC.canonical_frequency_collision(collision)
end

@testset "generated mixed canonical frequency assembly" begin
    collision = generated_mixed_canonical_collision()
    mixed_parameter = KC.ParameterMonomial(:g) * KC.ParameterMonomial(:γ)
    @test parameters(collision) == mixed_parameter
    @test target_family(collision) === canonical_mixed_ϕ

    # This generated mixed gγ sector is already unshifted at the kinetic boundary, matching the
    # previously certified mixed-sector oracle. The generic assembler still preserves finite
    # shifts whenever they are present; the genuinely shifted γ² sectors remain the acceptance
    # workload for the later Trotter layer.
    @test isempty(KC.shifted_frequency_terms(collision))

    expressions = KC.canonical_frequency_expressions(collision)
    @test !isempty(expressions)
    @test all(!isempty, values(expressions))

    census = NamedTuple[]
    for (sector, expression) in expressions
        nfrequencies = length(KC.momentum_basis(sector)) - 1
        decay = Tuple(KC.frequency_decay_lower_bound(expression, i) for i in 1:nfrequencies)
        push!(
            census,
            (;
                statistical=Tuple(
                    Tuple(KC.momentum(atom).coefficients) for
                    atom in KC.statistical_monomial(sector)
                ),
                nterms=length(KC.causal_frequency_terms(expression)),
                decay,
            ),
        )
    end
    sort!(census; by=repr)
    @info "generated mixed canonical frequency census" census

    # At least one generated sector must demonstrate why contour safety is an expression-level
    # property: the complete grouped expression is O(ω^-2) or better in some loop frequency.
    @test any(record -> any(>=(2), record.decay), census)
end
