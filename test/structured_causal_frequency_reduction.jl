using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields structured_causal_ϕ::Boson

function structured_causal_energies()
    basis = KC.MomentumBasis(2)
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    return EnergyForm(DispersionAtom(structured_causal_ϕ, k)),
    EnergyForm(DispersionAtom(structured_causal_ϕ, q))
end

function structured_causal_denominator(coefficient, energy, infinitesimal)
    return KC.CausalFrequencyDenominator([coefficient // 1], energy, infinitesimal // 1)
end

function structured_frequency_independent_denominator(energy, infinitesimal)
    return KC.CausalFrequencyDenominator([0 // 1], energy, infinitesimal // 1)
end

@testset "coincident causal pole witness" begin
    Eₖ, E_q = structured_causal_energies()

    upper = structured_causal_denominator(1, -E_q, -1)
    proportional_upper = structured_causal_denominator(2, -2 * E_q, -2)
    spectator = structured_causal_denominator(-1, Eₖ, -1)
    term = KC.CausalFrequencyTerm(3 // 2, [upper, proportional_upper, spectator])

    witness = @inferred KC.coincident_causal_pole_witness(term, 1)
    @test KC.has_coincident_causal_pole(witness)
    @test witness.frequency_index == 1
    @test witness.pivot_denominator != witness.coincident_denominator
    @test_throws ArgumentError KC.integrate_causal_frequency(term, 1)
end

@testset "opposite-side poles are not same-prescription coincidence" begin
    _, E_q = structured_causal_energies()
    retarded = structured_causal_denominator(1, -E_q, 1)
    advanced = structured_causal_denominator(1, -E_q, -1)
    term = KC.CausalFrequencyTerm(1 // 1, [retarded, advanced])

    witness = @inferred KC.coincident_causal_pole_witness(term, 1)
    @test !KC.has_coincident_causal_pole(witness)
end

@testset "same-side distinct poles retain simple residue path" begin
    Eₖ, E_q = structured_causal_energies()
    first = structured_causal_denominator(1, -E_q, -1)
    second = structured_causal_denominator(1, -Eₖ, -1)
    term = KC.CausalFrequencyTerm(1 // 1, [first, second])

    witness = @inferred KC.coincident_causal_pole_witness(term, 1)
    @test !KC.has_coincident_causal_pole(witness)
    @test length(@inferred KC.integrate_causal_frequency(term, 1)) == 2
end

@testset "zero-energy causal denominators remain exceptional" begin
    _, E_q = structured_causal_energies()
    zero_energy = zero(E_q)

    # Once all loop-frequency coefficients vanish, E=0 leaves either 1/(±i0) or PV(1/0).
    # Both are singular strict-quasiparticle objects and must be preserved rather than lowered.
    for infinitesimal in (-1, 0, 1)
        singular = structured_frequency_independent_denominator(zero_energy, infinitesimal)
        term = KC.CausalFrequencyTerm(1 // 1, [singular])
        witness = @inferred KC.zero_energy_causal_denominator_witness(term)
        @test KC.has_zero_energy_causal_denominator(witness)
        @test witness.denominator_index == 1
        @test_throws ArgumentError KC.lower_causal_frequency_term(term)
    end

    regular = structured_frequency_independent_denominator(E_q, 1)
    regular_term = KC.CausalFrequencyTerm(1 // 1, [regular])
    witness = @inferred KC.zero_energy_causal_denominator_witness(regular_term)
    @test !KC.has_zero_energy_causal_denominator(witness)
    @test !isempty(@inferred KC.lower_causal_frequency_term(regular_term))
end
