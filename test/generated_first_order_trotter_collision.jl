using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields canonical_trotter_loss_ϕ::Boson

function generated_first_order_canonical_trotter_collision()
    c = canonical_trotter_loss_ϕ[Classical]
    q = canonical_trotter_loss_ϕ[Quantum]
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
        L, Val(1), Val(3); simplify=true, _set_reg_to_zero=false
    )
    GF = @inferred fourier_transform(G)
    ΣF = @inferred SelfEnergy(GF)
    ΣW = @inferred wigner_transform(ΣF; gradient_order=Val(0))
    KΣ = @inferred kinetic_expression(ΣW)
    off_shell = @inferred off_shell_collision_expression(KΣ)
    spectral = @inferred spectral_dispersive_collision(off_shell)
    return @inferred KC.canonical_frequency_collision(spectral)
end

canonical_trotter_momentum_signature(momentum) = Tuple(momentum.coefficients)

function canonical_trotter_external_signature(carrier)
    basis = KC.momentum_basis(carrier)
    external = external_wigner_momentum(carrier)
    for (basis_index, variable) in enumerate(basis.variables)
        isequal(variable, external) || continue
        return canonical_trotter_momentum_signature(KC.basis_momentum(basis, basis_index))
    end
    return error("external momentum is absent from generated basis")
end

function canonical_trotter_statistical_signature(monomial)
    signature = [
        canonical_trotter_momentum_signature(KC.momentum(atom)) for atom in monomial
    ]
    sort!(signature)
    return Tuple(signature)
end

@testset "generated first-order loss through canonical Trotter layer" begin
    collision = generated_first_order_canonical_trotter_collision()
    @test parameters(collision) == KC.ParameterMonomial(:γ)
    @test target_family(collision) === canonical_trotter_loss_ϕ

    reduction = @inferred KC.reduce_canonical_trotter_frequencies(collision)
    @test isempty(KC.canonical_trotter_boundary_expressions(reduction))
    @test isempty(KC.blocked_trotter_contributions(reduction))
    @test isempty(KC.unresolved_trotter_states(reduction))

    grouped = Dict{Tuple,KC.ComplexRationals}()
    external_signatures = Set{Tuple}()
    constants = KC.canonical_trotter_constants(reduction)
    @test length(constants) == 4

    for (sector, coefficient) in constants
        push!(external_signatures, canonical_trotter_external_signature(sector))
        signature = canonical_trotter_statistical_signature(KC.statistical_monomial(sector))
        grouped[signature] =
            get(grouped, signature, zero(KC.ComplexRationals)) + coefficient
    end

    @test length(external_signatures) == 1
    external_signature = only(external_signatures)
    one_factor = Tuple[key for key in keys(grouped) if length(key) == 1]
    @test length(one_factor) == 2
    internal_signature = only(filter(key -> only(key) != external_signature, one_factor))
    external_monomial = (external_signature,)
    mixed_monomial = Tuple(sort([external_signature, only(internal_signature)]))

    # Independent pair-loss oracle:
    #
    # I_package^γ / γ = 2(F_k - 1)(1 - F_q)
    #                 = -2 + 2F_q + 2F_k - 2F_kF_q.
    @test get(grouped, (), zero(KC.ComplexRationals)) == -2 // 1
    @test get(grouped, internal_signature, zero(KC.ComplexRationals)) == 2 // 1
    @test get(grouped, external_monomial, zero(KC.ComplexRationals)) == 2 // 1
    @test get(grouped, mixed_monomial, zero(KC.ComplexRationals)) == -2 // 1
    @test length(grouped) == 4
end
