using KeldyshContraction, Test
import KeldyshContraction as KC

function partial_frequency_line(
    family::FieldFamily{S}, routed::KC.LinearMomentum
) where {S<:KC.Statistics}
    return KineticLine{S}(
        family,
        KC.Bulk(1),
        KC.Bulk(2),
        Int8(0),
        routed,
        KineticSpectral,
        NoStatisticalWeight,
    )
end

function partial_frequency_term(
    basis::KC.MomentumBasis,
    line_kinds::Vector{Pair{KineticLine{S},SpectralDispersiveKind}},
    ::Val{E},
) where {S<:KC.Statistics,E}
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

function partial_reduction_recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        partial_reduction_recursively_concrete(eltype(T), seen) || return false
    elseif T <: AbstractDict
        partial_reduction_recursively_concrete(keytype(T), seen) || return false
        partial_reduction_recursively_concrete(valtype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        partial_reduction_recursively_concrete(FT, seen) || return false
    end
    return true
end

@qfields partial_ϕ::Boson partial_χ::Boson

@testset "partial spectral elimination exposes causal residual subspace" begin
    basis = KC.MomentumBasis(3)
    k = KC.basis_momentum(basis, 1)
    q₁ = KC.basis_momentum(basis, 2)
    q₂ = KC.basis_momentum(basis, 3)
    q₃ = q₁ + q₂ - k

    line₁ = partial_frequency_line(partial_ϕ, q₁)
    line₂ = partial_frequency_line(partial_ϕ, q₂)
    line₃ = partial_frequency_line(partial_ϕ, q₃)
    term = partial_frequency_term(
        basis,
        [
            line₁ => CollisionSpectral,
            line₂ => CollisionDispersive,
            line₃ => CollisionDispersive,
        ],
        Val(3),
    )

    @test spectral_frequency_rank(term) == 1
    @test loop_frequency_count(term) == 2
    reduction = @inferred KC.partial_spectral_frequency_reduction(term, partial_ϕ)
    @test reduction isa KC.PartialSpectralFrequencyReduction{Boson}
    @test reduction.pivot_loop_basis_indices == [2]
    @test reduction.residual_loop_basis_indices == [3]
    @test length(reduction.pivot_spectral_lines) == 1
    @test reduction.pivot_residual_coefficients == zeros(KC.EnergyCoefficient, 1, 1)
    @test reduction.factor == 1 // 1
    @test isempty(reduction.support.shells)
    @test isempty(reduction.support.principal_values)
    @test length(reduction.causal_terms) == 4
    @test partial_reduction_recursively_concrete(typeof(reduction))

    εq₁ = EnergyForm(DispersionAtom(partial_ϕ, q₁))
    @test reduction.pivot_frequency_offsets == [εq₁]

    reduced = @inferred KC.reduce_partial_spectral_frequency(reduction)
    εk = EnergyForm(DispersionAtom(partial_ϕ, k))
    εq₂ = EnergyForm(DispersionAtom(partial_ϕ, q₂))
    εq₃ = EnergyForm(DispersionAtom(partial_ϕ, q₃))
    shell, shell_factor = energy_shell(εk + εq₃ - εq₁ - εq₂)
    @test shell_factor == 1 // 1
    expected_support = FrequencySupport([shell], PrincipalValueSupport{Boson}[])
    @test reduced == Dict(expected_support => convert(KC.ComplexRationals, 1 // 4))
end

@testset "general reducer preserves full-rank fast path exactly" begin
    basis = KC.MomentumBasis(3)
    k = KC.basis_momentum(basis, 1)
    q₁ = KC.basis_momentum(basis, 2)
    q₂ = KC.basis_momentum(basis, 3)
    q₃ = q₁ + q₂ - k

    line₁ = partial_frequency_line(partial_ϕ, q₁)
    line₂ = partial_frequency_line(partial_ϕ, q₂)
    line₃ = partial_frequency_line(partial_χ, q₃)
    term = partial_frequency_term(
        basis,
        [
            line₁ => CollisionSpectral,
            line₂ => CollisionSpectral,
            line₃ => CollisionDispersive,
        ],
        Val(3),
    )

    fast = @inferred full_rank_frequency_reduction(term, partial_ϕ)
    general = @inferred KC.general_frequency_reduction(term, partial_ϕ)
    expected = Dict(
        frequency_support(fast) => convert(KC.ComplexRationals, frequency_factor(fast))
    )
    @test general == expected
end

@testset "dependent shell retains independent residual support" begin
    basis = KC.MomentumBasis(2)
    q = KC.basis_momentum(basis, 2)

    line_ϕ = partial_frequency_line(partial_ϕ, q)
    line_χ = partial_frequency_line(partial_χ, q)
    term = partial_frequency_term(
        basis,
        [
            line_ϕ => CollisionSpectral,
            line_ϕ => CollisionSpectral,
            line_χ => CollisionSpectral,
        ],
        Val(3),
    )

    analysis = @inferred KC.analyze_spectral_dependencies(term, partial_ϕ)
    @test KC.has_dependent_shell_support(analysis)
    @test KC.constraint_rank(analysis) == 2
    @test KC.constraint_count(analysis) == 3

    reduction = @inferred KC.dependent_spectral_frequency_reduction(term, partial_ϕ)
    @test partial_reduction_recursively_concrete(typeof(reduction))
    @test reduction.factor == 1 // 1
    @test length(reduction.support.shells) == 1
    @test isempty(reduction.support.principal_values)

    reduced = @inferred KC.reduce_partial_spectral_frequency(reduction)
    εϕ = EnergyForm(DispersionAtom(partial_ϕ, q))
    εχ = EnergyForm(DispersionAtom(partial_χ, q))
    shell, shell_factor = energy_shell(εϕ - εχ)
    @test shell_factor == 1 // 1
    expected_support = FrequencySupport([shell], PrincipalValueSupport{Boson}[])
    @test reduced == Dict(expected_support => one(KC.ComplexRationals))
end
