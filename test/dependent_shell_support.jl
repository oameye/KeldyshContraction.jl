using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields dependent_shell_ϕ::Boson

zero_energy() = KC.EnergyForm{Boson}(3)
constraint(coefficients) = KC.AffineSpectralConstraint(coefficients, zero_energy())

@testset "dependent shell support" begin
    @testset "rank-one repeated shell" begin
        reduction = @inferred KC.analyze_spectral_dependencies([
            constraint([1 // 1, 0 // 1]), constraint([-2 // 1, 0 // 1])
        ])
        support = KC.dependent_shell_support(reduction)
        @test KC.has_dependent_shell_support(reduction)
        @test KC.constraint_rank(reduction) == 1
        @test KC.constraint_count(reduction) == 2
        @test KC.repeated_shell_multiplicity(reduction) == 2
        @test KC.constraint_scale_factor(reduction) == 1 // 2
        @test support.dependency_rows == [[1 // 1]]
    end

    @testset "row-order and scale invariance" begin
        a = KC.analyze_spectral_dependencies([
            constraint([1 // 1, 0 // 1]), constraint([-2 // 1, 0 // 1])
        ])
        b = KC.analyze_spectral_dependencies([
            constraint([6 // 1, 0 // 1]), constraint([-3 // 1, 0 // 1])
        ])
        @test KC.dependent_shell_support(a) == KC.dependent_shell_support(b)
        @test KC.constraint_scale_factor(a) == 1 // 2
        @test KC.constraint_scale_factor(b) == 1 // 18
    end

    @testset "higher multiplicity" begin
        reduction = @inferred KC.analyze_spectral_dependencies([
            constraint([1 // 1, 0 // 1]),
            constraint([1 // 1, 0 // 1]),
            constraint([1 // 1, 0 // 1]),
        ])
        support = KC.dependent_shell_support(reduction)
        @test KC.constraint_rank(reduction) == 1
        @test KC.constraint_count(reduction) == 3
        @test KC.repeated_shell_multiplicity(reduction) == 3
        @test support.dependency_rows == [[1 // 1], [1 // 1]]
    end

    @testset "higher-rank dependency geometry" begin
        reduction = @inferred KC.analyze_spectral_dependencies([
            constraint([1 // 1, 0 // 1]),
            constraint([0 // 1, 1 // 1]),
            constraint([1 // 1, 1 // 1]),
        ])
        support = KC.dependent_shell_support(reduction)
        @test KC.has_dependent_shell_support(reduction)
        @test KC.constraint_rank(reduction) == 2
        @test KC.constraint_count(reduction) == 3
        @test KC.repeated_shell_multiplicity(reduction) == 0
        @test support.dependency_rows == [[1 // 1, 1 // 1]]
    end

    @testset "same routing with different energy is regular" begin
        p = KC.LinearMomentum([1, 0, 0])
        q = KC.LinearMomentum([0, 1, 0])
        εp = KC.EnergyForm(KC.DispersionAtom(dependent_shell_ϕ, p))
        εq = KC.EnergyForm(KC.DispersionAtom(dependent_shell_ϕ, q))
        first_constraint = KC.AffineSpectralConstraint([1 // 1], εp)
        second_constraint = KC.AffineSpectralConstraint([1 // 1], εq)
        reduction = @inferred KC.analyze_spectral_dependencies([
            first_constraint, second_constraint
        ])
        @test !KC.has_dependent_shell_support(reduction)
        @test KC.constraint_rank(reduction) == 2
        @test KC.constraint_count(reduction) == 2
        @test isempty(KC.dependent_shell_support(reduction).dependency_rows)
    end

    @testset "regular independent system" begin
        reduction = @inferred KC.analyze_spectral_dependencies([
            constraint([1 // 1, 0 // 1]), constraint([0 // 1, 1 // 1])
        ])
        support = KC.dependent_shell_support(reduction)
        @test !KC.has_dependent_shell_support(reduction)
        @test KC.constraint_rank(reduction) == 2
        @test KC.constraint_count(reduction) == 2
        @test isempty(support.dependency_rows)
        @test KC.constraint_scale_factor(reduction) == 1
    end

    @testset "explicit delta zero" begin
        reduction = @inferred KC.analyze_spectral_dependencies([
            constraint([0 // 1, 0 // 1])
        ])
        support = KC.dependent_shell_support(reduction)
        @test KC.has_dependent_shell_support(reduction)
        @test KC.constraint_rank(reduction) == 0
        @test KC.constraint_count(reduction) == 1
        @test support.dependency_rows == [KC.EnergyCoefficient[]]
    end
end
