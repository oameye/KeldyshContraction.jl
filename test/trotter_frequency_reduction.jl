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

function trotter_state(
    term::SpectralDispersiveTerm{S,E1,E2};
    coefficient::KC.ComplexRationals=one(KC.ComplexRationals),
) where {S<:KC.Statistics,E1,E2}
    shifted = KC.ShiftedFrequencyCollisionTerm{KC.ComplexRationals,S,E1,E2}(
        term, coefficient, KC.StatisticalMonomial{S}(), KC.ParameterMonomial(:λ)
    )
    return KC.TrotterFrequencyState(shifted)
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
    state = @inferred trotter_state(term)

    witness = @inferred KC.trotter_frequency_witness(state)
    @test KC.has_isolated_trotter_frequency(witness)
    @test witness.frequency_index == 1
    @test witness.loop_basis_index == 2
    @test witness.loop_coefficient == 1
    @test regularisation_shift(KC.trotter_frequency_line(state, witness)) == 1
    @test @inferred(KC.trotter_equal_time_factor(state, witness)) == -(1 // 2) * im
end

@testset "one-sided retarded and advanced equal-time limits" begin
    basis = KC.MomentumBasis(2)
    q = KC.basis_momentum(basis, 2)

    function equal_time_components(shift)
        line = trotter_line(trotter_ϕ, q; shift)
        spectral = trotter_state(trotter_term(basis, [line => CollisionSpectral], Val(1)))
        dispersive = trotter_state(
            trotter_term(basis, [line => CollisionDispersive], Val(1))
        )
        spectral_witness = KC.trotter_frequency_witness(spectral)
        dispersive_witness = KC.trotter_frequency_witness(dispersive)
        A = KC.trotter_equal_time_factor(spectral, spectral_witness)
        D = KC.trotter_equal_time_factor(dispersive, dispersive_witness)
        return A, D
    end

    A_plus, D_plus = equal_time_components(1)
    @test D_plus - (1 // 2) * im * A_plus == -im
    @test iszero(D_plus + (1 // 2) * im * A_plus)

    A_minus, D_minus = equal_time_components(-1)
    @test iszero(D_minus - (1 // 2) * im * A_minus)
    @test D_minus + (1 // 2) * im * A_minus == im
end

@testset "spectral equal-time factor and routed Jacobian" begin
    basis = KC.MomentumBasis(3)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    shifted = trotter_line(trotter_ϕ, -2q + r; shift=-2)
    spectator = trotter_line(trotter_χ, r)

    spectral = trotter_state(
        trotter_term(
            basis, [shifted => CollisionSpectral, spectator => CollisionSpectral], Val(2)
        ),
    )
    spectral_witness = @inferred KC.trotter_frequency_witness(spectral)
    @test spectral_witness.loop_coefficient == -2
    @test @inferred(KC.trotter_equal_time_factor(spectral, spectral_witness)) == 1 // 2

    dispersive = trotter_state(
        trotter_term(
            basis, [shifted => CollisionDispersive, spectator => CollisionSpectral], Val(2)
        ),
    )
    dispersive_witness = @inferred KC.trotter_frequency_witness(dispersive)
    @test @inferred(KC.trotter_equal_time_factor(dispersive, dispersive_witness)) ==
        (1 // 4) * im
end

@testset "local Trotter rule refuses coupled or underconstrained frequencies" begin
    basis = KC.MomentumBasis(3)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)

    shifted = trotter_line(trotter_ϕ, q; shift=1)
    shared = trotter_line(trotter_χ, q + r)
    coupled = trotter_state(
        trotter_term(
            basis, [shifted => CollisionDispersive, shared => CollisionSpectral], Val(2)
        ),
    )
    @test !KC.has_isolated_trotter_frequency(
        @inferred KC.trotter_frequency_witness(coupled)
    )

    doubly_isolated = trotter_line(trotter_ϕ, q + r; shift=1)
    underconstrained = trotter_state(
        trotter_term(basis, [doubly_isolated => CollisionDispersive], Val(1))
    )
    @test !KC.has_isolated_trotter_frequency(
        @inferred KC.trotter_frequency_witness(underconstrained)
    )
end

@testset "sequential equal-time integration preserves spatial loop basis" begin
    basis = KC.MomentumBasis(3)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    shifted_q = trotter_line(trotter_ϕ, q; shift=1)
    shifted_r = trotter_line(trotter_χ, r; shift=-1)
    state = trotter_state(
        trotter_term(
            basis,
            [shifted_q => CollisionDispersive, shifted_r => CollisionDispersive],
            Val(2),
        ),
    )

    first_witness = @inferred KC.trotter_frequency_witness(state)
    first_state = @inferred KC.eliminate_trotter_frequency(state, first_witness)
    @test length(KC.active_trotter_lines(first_state)) == 1
    @test length(KC.active_trotter_frequencies(first_state)) == 1
    @test isequal(KC.momentum_basis(first_state), basis)

    second_witness = @inferred KC.trotter_frequency_witness(first_state)
    final_state = @inferred KC.eliminate_trotter_frequency(first_state, second_witness)
    @test isempty(KC.active_trotter_lines(final_state))
    @test isempty(KC.active_trotter_frequencies(final_state))
    @test KC.source_coefficient(final_state) == 1 // 4
    @test KC.trotter_frequency_complete(final_state)
    @test isequal(KC.momentum_basis(final_state), basis)
    @test length(KC.momentum_basis(final_state)) == 3
end

@testset "residual causal expression retains canonical frequency coordinates" begin
    basis = KC.MomentumBasis(3)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    shifted = trotter_line(trotter_ϕ, q; shift=1)
    spectator = trotter_line(trotter_χ, r)
    state = trotter_state(
        trotter_term(
            basis, [shifted => CollisionSpectral, spectator => CollisionSpectral], Val(2)
        ),
    )

    reduced = @inferred KC.eliminate_trotter_frequency(
        state, KC.trotter_frequency_witness(state)
    )
    @test KC.active_trotter_frequencies(reduced) == [2]
    @test !KC.has_unresolved_trotter_frequency(reduced)
    expression = @inferred KC.trotter_residual_causal_expression(reduced, trotter_ϕ)
    @test !isempty(expression)
    @test all(KC.causal_frequency_terms(expression)) do causal_term
        all(KC.causal_frequency_denominators(causal_term)) do denominator
            length(denominator.loop_coefficients) == 2 &&
                iszero(denominator.loop_coefficients[1])
        end
    end
end

@testset "statistically weighted shifted line is not locally universal" begin
    basis = KC.MomentumBasis(2)
    q = KC.basis_momentum(basis, 2)
    shifted = trotter_line(trotter_ϕ, q; shift=1, statistical=DistributionWeight)
    state = trotter_state(trotter_term(basis, [shifted => CollisionSpectral], Val(1)))
    witness = @inferred KC.trotter_frequency_witness(state)
    @test KC.has_isolated_trotter_frequency(witness)
    @test_throws ArgumentError KC.trotter_equal_time_factor(state, witness)
end

@qfields trotter_ψ::Fermion

@testset "equal-time spectral factor is statistics generic" begin
    basis = KC.MomentumBasis(2)
    q = KC.basis_momentum(basis, 2)
    shifted = trotter_line(trotter_ψ, q; shift=-1)
    state = trotter_state(trotter_term(basis, [shifted => CollisionSpectral], Val(1)))
    witness = @inferred KC.trotter_frequency_witness(state)
    @test KC.has_isolated_trotter_frequency(witness)
    @test @inferred(KC.trotter_equal_time_factor(state, witness)) == 1
end
