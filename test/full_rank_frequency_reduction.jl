using KeldyshContraction, Test
import KeldyshContraction as KC

function full_rank_recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        full_rank_recursively_concrete(eltype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        full_rank_recursively_concrete(FT, seen) || return false
    end
    return true
end

function frequency_line(
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

function synthetic_frequency_term(
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

@qfields frequency_ϕ::Boson frequency_χ::Boson

@testset "full-rank all-spectral sunset leaves one energy shell" begin
    basis = KC.MomentumBasis(3)
    k = KC.basis_momentum(basis, 1)
    q₁ = KC.basis_momentum(basis, 2)
    q₂ = KC.basis_momentum(basis, 3)
    q₃ = q₁ + q₂ - k

    line₁ = frequency_line(frequency_ϕ, q₁; statistical=DistributionWeight)
    line₂ = frequency_line(frequency_ϕ, q₂; statistical=DistributionWeight)
    line₃ = frequency_line(frequency_ϕ, q₃)
    term = synthetic_frequency_term(
        basis,
        [
            line₁ => CollisionSpectral,
            line₂ => CollisionSpectral,
            line₃ => CollisionSpectral,
        ],
        Val(3),
    )

    @test @inferred(loop_frequency_count(term)) == 2
    @test @inferred(spectral_frequency_rank(term)) == 2
    reduction = @inferred full_rank_frequency_reduction(term, frequency_ϕ)
    @test reduction isa FullRankFrequencyReduction{Boson}
    @test loop_frequency_basis_indices(reduction) == [2, 3]
    @test length(loop_frequency_energies(reduction)) == 2
    @test length(spectral_pivot_lines(reduction)) == 2
    @test frequency_factor(reduction) == 1 // 1
    @test full_rank_recursively_concrete(typeof(reduction))

    εk = EnergyForm(DispersionAtom(frequency_ϕ, k))
    εq₁ = EnergyForm(DispersionAtom(frequency_ϕ, q₁))
    εq₂ = EnergyForm(DispersionAtom(frequency_ϕ, q₂))
    εq₃ = EnergyForm(DispersionAtom(frequency_ϕ, q₃))
    expected_shell, expected_factor = energy_shell(εk + εq₃ - εq₁ - εq₂)
    support = frequency_support(reduction)
    @test support.shells == [expected_shell]
    @test isempty(support.principal_values)
    @test expected_factor == 1 // 1

    again = @inferred full_rank_frequency_reduction(term, frequency_ϕ)
    @test reduction == again
    @test isequal(reduction, again)
    @test hash(reduction) == hash(again)
end

@testset "same solver produces a principal-value energy denominator" begin
    basis = KC.MomentumBasis(3)
    k = KC.basis_momentum(basis, 1)
    q₁ = KC.basis_momentum(basis, 2)
    q₂ = KC.basis_momentum(basis, 3)
    q₃ = q₁ + q₂ - k

    line₁ = frequency_line(frequency_ϕ, q₁; statistical=DistributionWeight)
    line₂ = frequency_line(frequency_ϕ, q₂; statistical=DistributionWeight)
    line₃ = frequency_line(frequency_χ, q₃)
    term = synthetic_frequency_term(
        basis,
        [
            line₁ => CollisionSpectral,
            line₂ => CollisionSpectral,
            line₃ => CollisionDispersive,
        ],
        Val(3),
    )

    reduction = @inferred full_rank_frequency_reduction(term, frequency_ϕ)
    @test spectral_frequency_rank(term) == 2
    @test isempty(frequency_support(reduction).shells)
    @test length(frequency_support(reduction).principal_values) == 1

    εk = EnergyForm(DispersionAtom(frequency_ϕ, k))
    εq₁ = EnergyForm(DispersionAtom(frequency_ϕ, q₁))
    εq₂ = EnergyForm(DispersionAtom(frequency_ϕ, q₂))
    εq₃ = EnergyForm(DispersionAtom(frequency_χ, q₃))
    expected_pv, expected_factor = principal_value_support(-εk + εq₁ + εq₂ - εq₃)
    @test frequency_support(reduction).principal_values == [expected_pv]
    @test frequency_factor(reduction) == expected_factor
end

@testset "delta Jacobian is exact" begin
    basis = KC.MomentumBasis(3)
    q₁ = KC.basis_momentum(basis, 2)
    q₂ = KC.basis_momentum(basis, 3)
    line₁ = frequency_line(frequency_ϕ, 2 * q₁)
    line₂ = frequency_line(frequency_ϕ, 3 * q₂)
    term = synthetic_frequency_term(
        basis, [line₁ => CollisionSpectral, line₂ => CollisionSpectral], Val(2)
    )

    reduction = @inferred full_rank_frequency_reduction(term, frequency_ϕ)
    @test frequency_factor(reduction) == 1 // 6
    solved = loop_frequency_energies(reduction)
    @test solved[1] == (1 // 2) * EnergyForm(DispersionAtom(frequency_ϕ, 2 * q₁))
    @test solved[2] == (1 // 3) * EnergyForm(DispersionAtom(frequency_ϕ, 3 * q₂))
    @test isempty(frequency_support(reduction).shells)
    @test isempty(frequency_support(reduction).principal_values)
end

@testset "exceptional sectors are never silently reduced" begin
    basis = KC.MomentumBasis(3)
    q₁ = KC.basis_momentum(basis, 2)
    q₂ = KC.basis_momentum(basis, 3)
    line₁ = frequency_line(frequency_ϕ, q₁)
    line₂ = frequency_line(frequency_ϕ, q₂)
    line₃ = frequency_line(frequency_ϕ, q₁ + q₂)
    rank_deficient = synthetic_frequency_term(
        basis,
        [
            line₁ => CollisionSpectral,
            line₂ => CollisionDispersive,
            line₃ => CollisionDispersive,
        ],
        Val(3),
    )
    @test spectral_frequency_rank(rank_deficient) == 1
    @test_throws ArgumentError full_rank_frequency_reduction(rank_deficient, frequency_ϕ)

    one_loop_basis = KC.MomentumBasis(2)
    q = KC.basis_momentum(one_loop_basis, 2)
    spectral = frequency_line(frequency_ϕ, q)
    dispersive = frequency_line(frequency_ϕ, q)
    matched_ad = synthetic_frequency_term(
        one_loop_basis,
        [spectral => CollisionSpectral, dispersive => CollisionDispersive],
        Val(2),
    )
    @test spectral_frequency_rank(matched_ad) == 1
    @test_throws ArgumentError full_rank_frequency_reduction(matched_ad, frequency_ϕ)

    shifted = frequency_line(frequency_ϕ, q; shift=1)
    shifted_term = synthetic_frequency_term(
        one_loop_basis, [shifted => CollisionSpectral], Val(1)
    )
    @test_throws ArgumentError full_rank_frequency_reduction(shifted_term, frequency_ϕ)
end

@qfields frequency_ψ::Fermion

@testset "full-rank frequency reduction is statistics-generic" begin
    basis = KC.MomentumBasis(2)
    q = KC.basis_momentum(basis, 2)
    line = frequency_line(frequency_ψ, q; statistical=DistributionWeight)
    term = synthetic_frequency_term(basis, [line => CollisionSpectral], Val(1))
    reduction = @inferred full_rank_frequency_reduction(term, frequency_ψ)

    @test reduction isa FullRankFrequencyReduction{Fermion}
    @test spectral_frequency_rank(term) == 1
    @test frequency_factor(reduction) == 1 // 1
    @test loop_frequency_energies(reduction) == [EnergyForm(DispersionAtom(frequency_ψ, q))]
    @test full_rank_recursively_concrete(typeof(reduction))
end
