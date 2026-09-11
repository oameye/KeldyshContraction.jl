import KeldyshContraction as KC

@qfields benchmark_frequency_ϕ::Boson

function benchmark_frequency_line(routed::KC.LinearMomentum)
    return KineticLine{Boson}(
        benchmark_frequency_ϕ,
        KC.Bulk(1),
        KC.Bulk(2),
        Int8(0),
        routed,
        KineticSpectral,
        NoStatisticalWeight,
    )
end

function benchmark_frequency_term(
    basis::KC.MomentumBasis,
    line_kinds::Vector{Pair{KineticLine{Boson},SpectralDispersiveKind}},
    ::Val{E},
) where {E}
    canonical = copy(line_kinds)
    sort!(canonical; by=first)
    lines = KineticLine{Boson}[first(pair) for pair in canonical]
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

function benchmark_frequency_reduction!(suite)
    basis = KC.MomentumBasis(3)
    k = KC.basis_momentum(basis, 1)
    q₁ = KC.basis_momentum(basis, 2)
    q₂ = KC.basis_momentum(basis, 3)
    q₃ = q₁ + q₂ - k

    line₁ = benchmark_frequency_line(q₁)
    line₂ = benchmark_frequency_line(q₂)
    line₃ = benchmark_frequency_line(q₃)

    full_rank = benchmark_frequency_term(
        basis,
        [
            line₁ => CollisionSpectral,
            line₂ => CollisionSpectral,
            line₃ => CollisionSpectral,
        ],
        Val(3),
    )
    partial_causal = benchmark_frequency_term(
        basis,
        [
            line₁ => CollisionSpectral,
            line₂ => CollisionDispersive,
            line₃ => CollisionDispersive,
        ],
        Val(3),
    )
    dependent_shell = benchmark_frequency_term(
        basis,
        [
            line₁ => CollisionSpectral,
            line₁ => CollisionSpectral,
            line₂ => CollisionSpectral,
        ],
        Val(3),
    )

    suite["Frequency reduction"]["full-rank fast path"] = @benchmarkable KC.general_frequency_reduction(
        $full_rank, $benchmark_frequency_ϕ
    ) seconds = 10
    suite["Frequency reduction"]["partial causal path"] = @benchmarkable KC.general_frequency_reduction(
        $partial_causal, $benchmark_frequency_ϕ
    ) seconds = 10
    suite["Frequency reduction"]["dependent affine path"] = @benchmarkable KC.reduce_frequency_term(
        $dependent_shell, $benchmark_frequency_ϕ
    ) seconds = 10
    return suite
end
