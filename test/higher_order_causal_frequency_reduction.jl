using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields higher_order_causal_ϕ::Boson

function higher_order_causal_energies()
    basis = KC.MomentumBasis(2)
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    return EnergyForm(DispersionAtom(higher_order_causal_ϕ, k)),
    EnergyForm(DispersionAtom(higher_order_causal_ϕ, q))
end

function higher_order_denominator(coefficient, energy, infinitesimal)
    return KC.CausalFrequencyDenominator([coefficient // 1], energy, infinitesimal // 1)
end

@testset "causal pole classes use full normalized affine denominator" begin
    Eₖ, E_q = higher_order_causal_energies()
    first = higher_order_denominator(2, -2 * E_q, -2)
    proportional = higher_order_denominator(-3, 3 * E_q, 3)
    distinct = higher_order_denominator(1, -Eₖ, -1)
    term = KC.CausalFrequencyTerm(1 // 1, [first, proportional, distinct])

    classes = @inferred KC.causal_pole_classes(term, 1)
    @test sort!(length.(classes)) == [1, 2]
end

@testset "simple pole path is unchanged" begin
    Eₖ, E_q = higher_order_causal_energies()
    upper = higher_order_denominator(1, -E_q, -1)
    lower = higher_order_denominator(1, -Eₖ, 1)
    term = KC.CausalFrequencyTerm(3 // 2, [upper, lower])

    @test @inferred(KC.integrate_causal_frequency_exact(term, 1)) ==
        @inferred(KC.integrate_causal_frequency(term, 1))
end

@testset "pure same-side double pole integrates to zero" begin
    _, E_q = higher_order_causal_energies()
    upper = higher_order_denominator(1, -E_q, -1)
    term = KC.CausalFrequencyTerm(1 // 1, [upper, upper])

    @test isempty(@inferred KC.integrate_causal_frequency_exact(term, 1))
end

@testset "double-pole derivative residue" begin
    Eₖ, E_q = higher_order_causal_energies()
    upper = higher_order_denominator(1, -E_q, -1)
    lower = higher_order_denominator(1, -Eₖ, 1)
    term = KC.CausalFrequencyTerm(1 // 1, [upper, upper, lower])

    result = @inferred KC.integrate_causal_frequency_exact(term, 1)
    @test length(result) == 1
    transformed = higher_order_denominator(0, E_q - Eₖ, 2)
    expected = KC.CausalFrequencyTerm(
        convert(KC.ComplexRationals, -im), [transformed, transformed]
    )
    @test only(result) == expected
end

@testset "triple-pole derivative residue" begin
    Eₖ, E_q = higher_order_causal_energies()
    upper = higher_order_denominator(1, -E_q, -1)
    lower = higher_order_denominator(1, -Eₖ, 1)
    term = KC.CausalFrequencyTerm(1 // 1, [upper, upper, upper, lower])

    result = @inferred KC.integrate_causal_frequency_exact(term, 1)
    @test length(result) == 1
    transformed = higher_order_denominator(0, E_q - Eₖ, 2)
    expected = KC.CausalFrequencyTerm(
        convert(KC.ComplexRationals, im), [transformed, transformed, transformed]
    )
    @test only(result) == expected
end

@testset "non-unit proportional coefficients preserve exact Jacobian sign" begin
    Eₖ, E_q = higher_order_causal_energies()
    first = higher_order_denominator(2, -2 * E_q, -2)
    second = higher_order_denominator(-3, 3 * E_q, 3)
    lower = higher_order_denominator(5, -5 * Eₖ, 5)
    term = KC.CausalFrequencyTerm(1 // 1, [first, second, lower])

    result = @inferred KC.integrate_causal_frequency_exact(term, 1)
    @test length(result) == 1
    transformed = higher_order_denominator(0, 5 * (E_q - Eₖ), 10)
    expected = KC.CausalFrequencyTerm(
        convert(KC.ComplexRationals, (5 // 6) * im), [transformed, transformed]
    )
    @test only(result) == expected
end

@testset "repeated differentiation combines identical causal monomials" begin
    Eₖ, E_q = higher_order_causal_energies()
    upper = higher_order_denominator(1, -E_q, -1)
    lower = higher_order_denominator(1, -Eₖ, 1)
    term = KC.CausalFrequencyTerm(1 // 1, [upper, upper, lower, lower])

    result = @inferred KC.integrate_causal_frequency_exact(term, 1)
    @test length(result) == 1
    transformed = higher_order_denominator(0, E_q - Eₖ, 2)
    expected = KC.CausalFrequencyTerm(
        convert(KC.ComplexRationals, -2im), [transformed, transformed, transformed]
    )
    @test only(result) == expected
end

@testset "denominator permutation is canonical" begin
    Eₖ, E_q = higher_order_causal_energies()
    upper = higher_order_denominator(1, -E_q, -1)
    lower = higher_order_denominator(1, -Eₖ, 1)
    first = KC.CausalFrequencyTerm(2 // 3, [upper, upper, lower])
    second = KC.CausalFrequencyTerm(2 // 3, [lower, upper, upper])

    @test @inferred(KC.integrate_causal_frequency_exact(first, 1)) ==
        @inferred(KC.integrate_causal_frequency_exact(second, 1))
end

@testset "opposite prescriptions at one real pole remain distinct" begin
    _, E_q = higher_order_causal_energies()
    lower = higher_order_denominator(1, -E_q, 1)
    upper = higher_order_denominator(1, -E_q, -1)
    term = KC.CausalFrequencyTerm(1 // 1, [lower, upper])

    classes = @inferred KC.causal_pole_classes(term, 1)
    @test length(classes) == 2
    @test all(pole_class -> length(pole_class) == 1, classes)

    result = @inferred KC.integrate_causal_frequency_exact(term, 1)
    @test length(result) == 1
    @test KC.has_zero_energy_causal_denominator(
        KC.zero_energy_causal_denominator_witness(only(result))
    )
end
