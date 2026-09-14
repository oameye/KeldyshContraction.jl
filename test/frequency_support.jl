using KeldyshContraction, Test
import KeldyshContraction as KC

function frequency_recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        frequency_recursively_concrete(eltype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        frequency_recursively_concrete(FT, seen) || return false
    end
    return true
end

@qfields energy_ϕ::Boson energy_χ::Boson

@testset "canonical symbolic energy forms" begin
    p = KC.LinearMomentum([1, 0, 0])
    q = KC.LinearMomentum([0, 1, 0])
    r = KC.LinearMomentum([0, 0, 1])
    εp = DispersionAtom(energy_ϕ, p)
    εq = DispersionAtom(energy_ϕ, q)
    εrχ = DispersionAtom(energy_χ, r)

    form = @inferred EnergyForm(
        3, [εq => -1 // 1, εp => 1 // 1, εp => 2 // 1, εrχ => 0 // 1]
    )
    expected = EnergyForm(3, [εp => 3 // 1, εq => -1 // 1])
    @test form == expected
    @test energy_basis_size(form) == 3
    @test length(form) == 2
    @test collect(energy_terms(form)) == collect(energy_terms(expected))
    @test frequency_recursively_concrete(typeof(form))

    @test form - form == zero(form)
    @test iszero(form - form)
    @test 2 * EnergyForm(εp) - EnergyForm(εp) == EnergyForm(εp)

    wrong_dimension = DispersionAtom(energy_ϕ, KC.LinearMomentum([1, 0]))
    @test_throws DimensionMismatch EnergyForm(3, [wrong_dimension => 1 // 1])
    @test DispersionAtom(energy_ϕ, p) != DispersionAtom(energy_χ, p)
end

@testset "shell and principal-value rational-scale normal forms" begin
    p = KC.LinearMomentum([1, 0])
    q = KC.LinearMomentum([0, 1])
    εp = DispersionAtom(energy_ϕ, p)
    εq = DispersionAtom(energy_ϕ, q)
    mismatch = EnergyForm(εp) - EnergyForm(εq)

    shell, shell_factor = @inferred energy_shell(mismatch)
    reversed_shell, reversed_shell_factor = @inferred energy_shell(-mismatch)
    scaled_shell, scaled_shell_factor = @inferred energy_shell((3 // 2) * mismatch)
    @test shell == reversed_shell == scaled_shell
    @test isequal(shell, reversed_shell)
    @test hash(shell) == hash(reversed_shell)
    @test shell_factor == 1 // 1
    @test reversed_shell_factor == 1 // 1
    @test scaled_shell_factor == 2 // 3
    @test frequency_recursively_concrete(typeof(shell))

    pv, pv_factor = @inferred principal_value_support(mismatch)
    reversed_pv, reversed_pv_factor = @inferred principal_value_support(-mismatch)
    scaled_pv, scaled_pv_factor = @inferred principal_value_support((-3 // 2) * mismatch)
    @test pv == reversed_pv == scaled_pv
    @test pv_factor * mismatch == pv.energy
    @test reversed_pv_factor * (-mismatch) == reversed_pv.energy
    @test scaled_pv_factor * ((-3 // 2) * mismatch) == scaled_pv.energy
    @test frequency_recursively_concrete(typeof(pv))

    @test_throws ArgumentError energy_shell(zero(mismatch))
    @test_throws ArgumentError principal_value_support(zero(mismatch))
end

@testset "frequency support has deterministic factor ordering" begin
    p = KC.LinearMomentum([1, 0, 0])
    q = KC.LinearMomentum([0, 1, 0])
    r = KC.LinearMomentum([0, 0, 1])
    εp = EnergyForm(DispersionAtom(energy_ϕ, p))
    εq = EnergyForm(DispersionAtom(energy_ϕ, q))
    εr = EnergyForm(DispersionAtom(energy_ϕ, r))

    shell₁, _ = energy_shell(εp - εq)
    shell₂, _ = energy_shell(εp - εr)
    pv₁, _ = principal_value_support(εq - εr)
    pv₂, _ = principal_value_support(εp - εr)

    support_a = @inferred FrequencySupport([shell₂, shell₁], [pv₂, pv₁])
    support_b = @inferred FrequencySupport([shell₁, shell₂], [pv₁, pv₂])
    @test support_a == support_b
    @test isequal(support_a, support_b)
    @test hash(support_a) == hash(support_b)
    @test frequency_recursively_concrete(typeof(support_a))
end
