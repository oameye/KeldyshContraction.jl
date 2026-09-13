using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields kk_precedence_ϕ::Boson

function kk_precedence_line(
    routed::KC.LinearMomentum; statistical::StatisticalWeight=NoStatisticalWeight
)
    return KineticLine{Boson}(
        kk_precedence_ϕ,
        KC.Bulk(1),
        KC.Bulk(2),
        Int8(0),
        routed,
        KineticSpectral,
        statistical,
    )
end

function kk_precedence_term()
    basis = KC.MomentumBasis(3)
    q₁ = KC.basis_momentum(basis, 2)
    q₂ = KC.basis_momentum(basis, 3)

    spectral_q₁ = kk_precedence_line(q₁; statistical=DistributionWeight)
    dispersive_q₁ = kk_precedence_line(q₁)
    spectral_q₂_a = kk_precedence_line(q₂; statistical=DistributionWeight)
    spectral_q₂_b = kk_precedence_line(q₂)

    line_kinds = [
        spectral_q₁ => CollisionSpectral,
        dispersive_q₁ => CollisionDispersive,
        spectral_q₂_a => CollisionSpectral,
        spectral_q₂_b => CollisionSpectral,
    ]
    sort!(line_kinds; by=first)
    lines = KineticLine{Boson}[first(pair) for pair in line_kinds]
    kinds = SpectralDispersiveKind[last(pair) for pair in line_kinds]
    carrier = KineticTerm(
        KineticMonomial(lines, Val(4)),
        KC.FixedVector{0,Int}(Int[]),
        basis,
        basis[1],
        KC._kinematic_identity(),
    )
    return SpectralDispersiveTerm(
        carrier, KC.FixedVector{4,SpectralDispersiveKind}(kinds)
    )
end

@testset "Kramers-Kronig zero precedes dependent shell preservation" begin
    term = kk_precedence_term()
    dependency = @inferred KC.analyze_spectral_dependencies(term, kk_precedence_ϕ)
    @test KC.has_dependent_shell_support(dependency)

    classification = @inferred classify_exceptional_frequency(term)
    @test exceptional_frequency_kind(classification) === FrequencyKramersKronigZero

    # The isolated A(q₁)D(q₁) integral annihilates the complete term. The duplicated q₂ shell
    # is therefore never a physical pinch contribution and must not be preserved first.
    result = @inferred KC.reduce_frequency_term(term, kk_precedence_ϕ)
    @test isempty(result)
end
