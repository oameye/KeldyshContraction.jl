using KeldyshContraction, Test
import KeldyshContraction as KC

function exceptional_recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        exceptional_recursively_concrete(eltype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        exceptional_recursively_concrete(FT, seen) || return false
    end
    return true
end

function exceptional_line(
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

function exceptional_term(
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

@qfields exceptional_ϕ::Boson exceptional_χ::Boson

@testset "isolated spectral/dispersive pair proves Kramers-Kronig zero" begin
    basis = KC.MomentumBasis(2)
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    spectral = exceptional_line(exceptional_ϕ, q; statistical=DistributionWeight)
    dispersive = exceptional_line(exceptional_ϕ, q)
    external_only = exceptional_line(exceptional_χ, k)
    term = exceptional_term(
        basis,
        [
            spectral => CollisionSpectral,
            dispersive => CollisionDispersive,
            external_only => CollisionSpectral,
        ],
        Val(3),
    )

    classification = @inferred classify_exceptional_frequency(term)
    @test classification isa ExceptionalFrequencyClassification{Boson,3,0}
    @test exceptional_frequency_kind(classification) === FrequencyKramersKronigZero
    @test exceptional_frequency_term(classification) == term
    @test classification.spectral_rank == 1
    @test classification.loop_count == 1
    @test exceptional_frequency_loop_basis_index(classification) == 2

    spectral_line, dispersive_line = exceptional_frequency_witness_lines(classification)
    @test spectral_dispersive_kinds(term)[spectral_line] === CollisionSpectral
    @test spectral_dispersive_kinds(term)[dispersive_line] === CollisionDispersive
    @test KC.momentum(kinetic_lines(term.carrier)[spectral_line]) == q
    @test KC.momentum(kinetic_lines(term.carrier)[dispersive_line]) == q
    @test exceptional_recursively_concrete(typeof(classification))

    again = @inferred classify_exceptional_frequency(term)
    @test classification == again
    @test isequal(classification, again)
    @test hash(classification) == hash(again)
end

@testset "matched A-D singularity is handed to the exceptional rule" begin
    basis = KC.MomentumBasis(2)
    q = KC.basis_momentum(basis, 2)
    spectral = exceptional_line(exceptional_ϕ, q)
    dispersive = exceptional_line(exceptional_ϕ, q)
    term = exceptional_term(
        basis, [spectral => CollisionSpectral, dispersive => CollisionDispersive], Val(2)
    )

    @test spectral_frequency_rank(term) == loop_frequency_count(term) == 1
    @test_throws ArgumentError full_rank_frequency_reduction(term, exceptional_ϕ)
    classification = @inferred classify_exceptional_frequency(term)
    @test exceptional_frequency_kind(classification) === FrequencyKramersKronigZero
end

@testset "Kramers-Kronig proof refuses coupled and unmatched frequency structure" begin
    basis = KC.MomentumBasis(2)
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    spectral = exceptional_line(exceptional_ϕ, q; statistical=DistributionWeight)
    dispersive = exceptional_line(exceptional_ϕ, q)
    coupled = exceptional_line(exceptional_χ, q + k)
    coupled_term = exceptional_term(
        basis,
        [
            spectral => CollisionSpectral,
            dispersive => CollisionDispersive,
            coupled => CollisionSpectral,
        ],
        Val(3),
    )
    coupled_classification = @inferred classify_exceptional_frequency(coupled_term)
    @test exceptional_frequency_kind(coupled_classification) === FrequencyUnresolved
    @test exceptional_frequency_witness_lines(coupled_classification) == (0, 0)
    @test exceptional_frequency_loop_basis_index(coupled_classification) == 0

    unmatched = exceptional_line(exceptional_ϕ, q + k)
    unmatched_term = exceptional_term(
        basis, [spectral => CollisionSpectral, unmatched => CollisionDispersive], Val(2)
    )
    unmatched_classification = @inferred classify_exceptional_frequency(unmatched_term)
    @test exceptional_frequency_kind(unmatched_classification) === FrequencyUnresolved
end

@testset "equal-time shifts take precedence over causal zero classification" begin
    basis = KC.MomentumBasis(2)
    q = KC.basis_momentum(basis, 2)
    shifted_spectral = exceptional_line(
        exceptional_ϕ, q; shift=1, statistical=DistributionWeight
    )
    shifted_dispersive = exceptional_line(exceptional_ϕ, q; shift=1)
    term = exceptional_term(
        basis,
        [shifted_spectral => CollisionSpectral, shifted_dispersive => CollisionDispersive],
        Val(2),
    )

    classification = @inferred classify_exceptional_frequency(term)
    @test exceptional_frequency_kind(classification) === FrequencyTrotterRequired
    shifted_line, second = exceptional_frequency_witness_lines(classification)
    @test shifted_line > 0
    @test second == 0
    @test !iszero(regularisation_shift(kinetic_lines(term.carrier)[shifted_line]))
    @test exceptional_frequency_loop_basis_index(classification) == 0
end

@testset "rank deficiency alone never implies zero or principal value" begin
    basis = KC.MomentumBasis(3)
    q₁ = KC.basis_momentum(basis, 2)
    q₂ = KC.basis_momentum(basis, 3)
    spectral = exceptional_line(exceptional_ϕ, q₁)
    dispersive₁ = exceptional_line(exceptional_ϕ, q₂)
    dispersive₂ = exceptional_line(exceptional_χ, q₁ + q₂)
    term = exceptional_term(
        basis,
        [
            spectral => CollisionSpectral,
            dispersive₁ => CollisionDispersive,
            dispersive₂ => CollisionDispersive,
        ],
        Val(3),
    )

    @test spectral_frequency_rank(term) == 1
    @test loop_frequency_count(term) == 2
    classification = @inferred classify_exceptional_frequency(term)
    @test exceptional_frequency_kind(classification) === FrequencyUnresolved
    @test classification.spectral_rank == 1
    @test classification.loop_count == 2
end

@qfields exceptional_ψ::Fermion

@testset "exceptional classification is statistics-generic" begin
    basis = KC.MomentumBasis(2)
    q = KC.basis_momentum(basis, 2)
    spectral = exceptional_line(exceptional_ψ, q; statistical=DistributionWeight)
    dispersive = exceptional_line(exceptional_ψ, q)
    term = exceptional_term(
        basis, [spectral => CollisionSpectral, dispersive => CollisionDispersive], Val(2)
    )

    classification = @inferred classify_exceptional_frequency(term)
    @test classification isa ExceptionalFrequencyClassification{Fermion,2,0}
    @test exceptional_frequency_kind(classification) === FrequencyKramersKronigZero
    @test exceptional_recursively_concrete(typeof(classification))
end
