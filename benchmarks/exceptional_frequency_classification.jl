import KeldyshContraction as KC

@qfields benchmark_exceptional_ϕ::Boson

function benchmark_exceptional_frequency_classification!(suite)
    basis = KC.MomentumBasis(2)
    q = KC.basis_momentum(basis, 2)
    spectral = KineticLine{Boson}(
        benchmark_exceptional_ϕ,
        KC.Bulk(1),
        KC.Bulk(2),
        Int8(0),
        q,
        KineticSpectral,
        DistributionWeight,
    )
    dispersive = KineticLine{Boson}(
        benchmark_exceptional_ϕ,
        KC.Bulk(1),
        KC.Bulk(2),
        Int8(0),
        q,
        KineticSpectral,
        NoStatisticalWeight,
    )
    lines = KineticLine{Boson}[dispersive, spectral]
    sort!(lines)
    kinds = SpectralDispersiveKind[
        if statistical_weight(line) === DistributionWeight
            CollisionSpectral
        else
            CollisionDispersive
        end for line in lines
    ]
    carrier = KineticTerm(
        KineticMonomial(lines, Val(2)),
        KC.FixedVector{0,Int}(Int[]),
        basis,
        basis[1],
        KC._kinematic_identity(),
    )
    term = SpectralDispersiveTerm(carrier, KC.FixedVector{2,SpectralDispersiveKind}(kinds))

    suite["Frequency reduction"]["exceptional classification"] = @benchmarkable classify_exceptional_frequency(
        $term
    ) seconds = 10
    return suite
end
