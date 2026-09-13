using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields causal_ϕ::Boson

function causal_test_energies()
    basis = KC.MomentumBasis(3)
    momenta = ntuple(i -> KC.basis_momentum(basis, i), 3)
    return map(momentum -> EnergyForm(DispersionAtom(causal_ϕ, momentum)), momenta)
end

function causal_denominator(coefficient, energy, infinitesimal)
    return KC.CausalFrequencyDenominator([coefficient // 1], energy, infinitesimal // 1)
end

function add_supports!(out, values)
    for (support, coefficient) in values
        combined = get(out, support, zero(KC.ComplexRationals)) + coefficient
        if iszero(combined)
            delete!(out, support)
        else
            out[support] = combined
        end
    end
    return out
end

function integrate_and_lower!(out, term)
    residues = @inferred KC.integrate_causal_frequency(term, 1)
    for residue in residues
        add_supports!(out, @inferred(KC.lower_causal_frequency_term(residue)))
    end
    return out
end

@testset "exact causal frequency reduction" begin
    Eₐ, Eᵦ, Eₖ = causal_test_energies()

    @testset "same-half-plane causal product vanishes" begin
        term = KC.CausalFrequencyTerm(
            1 // 1, [causal_denominator(1, -Eₐ, 1), causal_denominator(1, -Eᵦ, 1)]
        )
        @test isempty(@inferred KC.integrate_causal_frequency(term, 1))
    end

    @testset "opposing-frequency dispersive convolution generates a shell" begin
        reduced = Dict{FrequencySupport{Boson},KC.ComplexRationals}()
        for η₁ in (-1, 1), η₂ in (-1, 1)
            term = KC.CausalFrequencyTerm(
                1 // 4,
                [causal_denominator(1, -Eₐ, η₁), causal_denominator(-1, Eₖ - Eᵦ, η₂)],
            )
            integrate_and_lower!(reduced, term)
        end

        shell, _ = energy_shell(Eₖ - Eₐ - Eᵦ)
        expected = FrequencySupport([shell], PrincipalValueSupport{Boson}[])
        @test reduced == Dict(expected => convert(KC.ComplexRationals, -1 // 4))
    end

    @testset "same-frequency dispersive convolution generates the opposite shell sign" begin
        reduced = Dict{FrequencySupport{Boson},KC.ComplexRationals}()
        for η₁ in (-1, 1), η₂ in (-1, 1)
            term = KC.CausalFrequencyTerm(
                1 // 4, [causal_denominator(1, -Eₐ, η₁), causal_denominator(1, -Eᵦ, η₂)]
            )
            integrate_and_lower!(reduced, term)
        end

        shell, _ = energy_shell(Eₐ - Eᵦ)
        expected = FrequencySupport([shell], PrincipalValueSupport{Boson}[])
        @test reduced == Dict(expected => convert(KC.ComplexRationals, 1 // 4))
    end

    @testset "non-unit frequency coefficients retain exact residue Jacobian" begin
        pivot = causal_denominator(2, -Eₐ, -1)
        other = causal_denominator(3, -Eᵦ, 1)
        term = KC.CausalFrequencyTerm(1 // 1, [pivot, other])
        residues = @inferred KC.integrate_causal_frequency(term, 1)

        @test length(residues) == 1
        residue = only(residues)
        @test residue.coefficient == convert(KC.ComplexRationals, (1 // 2) * im)
        @test length(residue.denominators) == 1
        denominator = only(residue.denominators)
        @test denominator.loop_coefficients == [0 // 1]
        @test denominator.energy == -Eᵦ + (3 // 2) * Eₐ
        @test denominator.infinitesimal == 5 // 2
    end
end
