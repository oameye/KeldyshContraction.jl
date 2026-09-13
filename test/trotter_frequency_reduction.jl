using KeldyshContraction, Test
import KeldyshContraction as KC

function trotter_line(
    family::FieldFamily{S},
    routed::KC.LinearMomentum;
    shift::Integer=0,
    statistical::StatisticalWeight=NoStatisticalWeight,
) where {S<:KC.Statistics}
    return KineticLine{S}(
        family,
        KC.Bulk(1),
        KC.Bulk(2),
        convert(Int8, shift),
        routed,
        KineticSpectral,
        statistical,
    )
end

function trotter_term(
    basis::KC.MomentumBasis,
    line_kinds::Vector{Pair{KineticLine{S},SpectralDispersiveKind}},
    ::Val{E},
) where {S<:KC.Statistics,E}
    length(line_kinds) == E || throw(ArgumentError("synthetic line count mismatch"))
    canonical = copy(line_kinds)
    sort!(canonical; by=first)
    lines = KineticLine{S}[first(pair) for pair in canonical]
    kinds = SpectralDispersiveKind[last(pair) for pair in canonical]
    carrier = KineticTerm(
        KineticMonomial(lines, Val(E)),
        KC.FixedVector{0,Int}(Int[]),
        basis,
        basis[1],
        KC._kinematic_identity(),
    )
    return SpectralDispersiveTerm(carrier, KC.FixedVector{E,SpectralDispersiveKind}(kinds))
end

@qfields trotter_ϕ::Boson trotter_χ::Boson

@testset "isolated shifted dispersive frequency" begin
    basis = KC.MomentumBasis(3)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    shifted = trotter_line(trotter_ϕ, q + r; shift=1)
    spectator = trotter_line(trotter_χ, r)
    term = trotter_term(
        basis, [shifted => CollisionDispersive, spectator => CollisionSpectral], Val(2)
    )

    witness = @inferred KC.trotter_frequency_witness(term)
    @test KC.has_isolated_trotter_frequency(witness)
    @test witness.loop_basis_index == 2
    @test witness.loop_coefficient == 1
    @test regularisation_shift(KC.trotter_frequency_line(term, witness)) == 1
    @test @inferred(KC.trotter_equal_time_factor(term, witness)) == -(1 // 2) * im
end

@testset "one-sided retarded and advanced equal-time limits" begin
    basis = KC.MomentumBasis(2)
    q = KC.basis_momentum(basis, 2)

    function equal_time_components(shift)
        line = trotter_line(trotter_ϕ, q; shift)
        spectral = trotter_term(basis, [line => CollisionSpectral], Val(1))
        dispersive = trotter_term(basis, [line => CollisionDispersive], Val(1))
        spectral_witness = @inferred KC.trotter_frequency_witness(spectral)
        dispersive_witness = @inferred KC.trotter_frequency_witness(dispersive)
        A = @inferred KC.trotter_equal_time_factor(spectral, spectral_witness)
        D = @inferred KC.trotter_equal_time_factor(dispersive, dispersive_witness)
        return A, D
    end

    A_plus, D_plus = equal_time_components(1)
    retarded_plus = D_plus - (1 // 2) * im * A_plus
    advanced_plus = D_plus + (1 // 2) * im * A_plus
    @test retarded_plus == -im
    @test iszero(advanced_plus)

    A_minus, D_minus = equal_time_components(-1)
    retarded_minus = D_minus - (1 // 2) * im * A_minus
    advanced_minus = D_minus + (1 // 2) * im * A_minus
    @test iszero(retarded_minus)
    @test advanced_minus == im
end

@testset "spectral equal-time factor and routed Jacobian" begin
    basis = KC.MomentumBasis(3)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    shifted = trotter_line(trotter_ϕ, -2q + r; shift=-2)
    spectator = trotter_line(trotter_χ, r)

    spectral_term = trotter_term(
        basis, [shifted => CollisionSpectral, spectator => CollisionSpectral], Val(2)
    )
    spectral_witness = @inferred KC.trotter_frequency_witness(spectral_term)
    @test KC.has_isolated_trotter_frequency(spectral_witness)
    @test spectral_witness.loop_coefficient == -2
    @test @inferred(KC.trotter_equal_time_factor(spectral_term, spectral_witness)) == 1 // 2

    dispersive_term = trotter_term(
        basis, [shifted => CollisionDispersive, spectator => CollisionSpectral], Val(2)
    )
    dispersive_witness = @inferred KC.trotter_frequency_witness(dispersive_term)
    @test @inferred(KC.trotter_equal_time_factor(dispersive_term, dispersive_witness)) ==
        (1 // 4) * im
end

@testset "local Trotter rule refuses coupled or underconstrained frequencies" begin
    basis = KC.MomentumBasis(3)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)

    shifted = trotter_line(trotter_ϕ, q; shift=1)
    shared = trotter_line(trotter_χ, q + r)
    coupled = trotter_term(
        basis, [shifted => CollisionDispersive, shared => CollisionSpectral], Val(2)
    )
    coupled_witness = @inferred KC.trotter_frequency_witness(coupled)
    @test !KC.has_isolated_trotter_frequency(coupled_witness)

    doubly_isolated = trotter_line(trotter_ϕ, q + r; shift=1)
    underconstrained = trotter_term(basis, [doubly_isolated => CollisionDispersive], Val(1))
    underconstrained_witness = @inferred KC.trotter_frequency_witness(underconstrained)
    @test !KC.has_isolated_trotter_frequency(underconstrained_witness)
end

@testset "sequential Trotter elimination preserves spatial loop basis" begin
    basis = KC.MomentumBasis(3)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    shifted_q = trotter_line(trotter_ϕ, q; shift=1)
    shifted_r = trotter_line(trotter_χ, r; shift=-1)
    term = trotter_term(
        basis, [shifted_q => CollisionDispersive, shifted_r => CollisionDispersive], Val(2)
    )

    state = @inferred KC.FrequencyIntegrationState(term)
    @test length(state.active_line_indices) == 2
    @test length(state.active_loop_basis_indices) == 2

    first_witness = @inferred KC.trotter_frequency_witness(state)
    first_state = @inferred KC.eliminate_trotter_frequency(state, first_witness)
    @test length(first_state.active_line_indices) == 1
    @test length(first_state.active_loop_basis_indices) == 1
    @test isequal(KC.momentum_basis(first_state.term), basis)

    second_witness = @inferred KC.trotter_frequency_witness(first_state)
    final_state = @inferred KC.eliminate_trotter_frequency(first_state, second_witness)
    @test isempty(final_state.active_line_indices)
    @test isempty(final_state.active_loop_basis_indices)
    @test final_state.factor == 1 // 4
    @test isequal(KC.momentum_basis(final_state.term), basis)
    @test length(KC.momentum_basis(final_state.term)) == 3
end

@testset "active dependency analysis ignores collapsed equal-time factor" begin
    basis = KC.MomentumBasis(3)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    shifted = trotter_line(trotter_ϕ, q; shift=1)
    shell₁ = trotter_line(trotter_χ, r)
    shell₂ = trotter_line(trotter_χ, r)
    term = trotter_term(
        basis,
        [
            shifted => CollisionDispersive,
            shell₁ => CollisionSpectral,
            shell₂ => CollisionSpectral,
        ],
        Val(3),
    )

    state = @inferred KC.FrequencyIntegrationState(term)
    witness = @inferred KC.trotter_frequency_witness(state)
    reduced = @inferred KC.eliminate_trotter_frequency(state, witness)
    analysis = @inferred KC.analyze_spectral_dependencies(reduced, trotter_ϕ)

    @test reduced.active_loop_basis_indices == [3]
    @test KC.has_dependent_shell_support(analysis)
    @test KC.constraint_rank(analysis) == 1
    @test KC.constraint_count(analysis) == 2
    @test KC.repeated_shell_multiplicity(analysis) == 2
    @test isequal(KC.momentum_basis(reduced.term), basis)
end

@testset "statistically weighted shifted line is not a universal equal-time constant" begin
    basis = KC.MomentumBasis(2)
    q = KC.basis_momentum(basis, 2)
    shifted = trotter_line(trotter_ϕ, q; shift=1, statistical=DistributionWeight)
    term = trotter_term(basis, [shifted => CollisionSpectral], Val(1))
    witness = @inferred KC.trotter_frequency_witness(term)
    @test KC.has_isolated_trotter_frequency(witness)
    @test_throws ArgumentError KC.trotter_equal_time_factor(term, witness)
end

@qfields trotter_ψ::Fermion

@testset "equal-time spectral factor is statistics generic" begin
    basis = KC.MomentumBasis(2)
    q = KC.basis_momentum(basis, 2)
    shifted = trotter_line(trotter_ψ, q; shift=-1)
    term = trotter_term(basis, [shifted => CollisionSpectral], Val(1))
    witness = @inferred KC.trotter_frequency_witness(term)
    @test KC.has_isolated_trotter_frequency(witness)
    @test @inferred(KC.trotter_equal_time_factor(term, witness)) == 1
end

@qfields trotter_loss_ϕ::Boson

function generated_first_order_trotter_collision()
    c = trotter_loss_ϕ[Classical]
    q = trotter_loss_ϕ[Quantum]
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

    L = InteractionLagrangian(loss, :γ)
    G = DressedPropagator(L, Val(1), Val(3); simplify=true, _set_reg_to_zero=false)
    GF = fourier_transform(G)
    ΣW = wigner_transform(SelfEnergy(GF); gradient_order=Val(0))
    collision = off_shell_collision_expression(kinetic_expression(ΣW))
    return spectral_dispersive_collision(collision)
end

@testset "generated first-order loss Trotter census" begin
    collision = generated_first_order_trotter_collision()
    @test parameters(collision) == KC.ParameterMonomial(:γ)

    shifted_terms = 0
    isolated_terms = 0
    dispersive_isolated = 0
    spectral_isolated = 0

    for expression in
        (collision_offset(collision), collision_distribution_coefficient(collision))
        for (term, _) in expression
            any(line -> !iszero(regularisation_shift(line)), kinetic_lines(term.carrier)) ||
                continue
            shifted_terms += 1
            witness = @inferred KC.trotter_frequency_witness(term)
            KC.has_isolated_trotter_frequency(witness) || continue
            isolated_terms += 1
            kind = spectral_dispersive_kinds(term)[witness.line_index]
            if kind === CollisionDispersive
                dispersive_isolated += 1
            else
                spectral_isolated += 1
            end
            @test statistical_weight(KC.trotter_frequency_line(term, witness)) ===
                NoStatisticalWeight
            @test @inferred(KC.trotter_equal_time_factor(term, witness)) isa
                KC.ComplexRationals
        end
    end

    @test shifted_terms == 4
    @test isolated_terms == 4
    @test dispersive_isolated == 2
    @test spectral_isolated == 2
end
