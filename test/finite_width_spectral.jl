using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields finite_width_ϕ::Boson

function finite_width_line(
    family::FieldFamily{S}, routed::KC.LinearMomentum
) where {S<:KC.Statistics}
    return KC.KineticLine{S}(
        family,
        KC.Bulk(1),
        KC.Bulk(2),
        Int8(0),
        routed,
        KC.KineticSpectral,
        KC.NoStatisticalWeight,
    )
end

function finite_width_term(
    basis::KC.MomentumBasis, lines::Vector{KC.KineticLine{S}}, ::Val{E}
) where {S<:KC.Statistics,E}
    length(lines) == E || throw(ArgumentError("synthetic line count mismatch"))
    carrier = KC.KineticTerm(
        KC.KineticMonomial(lines, Val(E)),
        KC.FixedVector{0,Int}(Int[]),
        basis,
        basis[1],
        KC._kinematic_identity(),
    )
    kinds = KC.FixedVector{E,KC.SpectralDispersiveKind}(fill(KC.CollisionSpectral, E))
    return KC.SpectralDispersiveTerm(carrier, kinds)
end

@testset "exact Lorentzian spectral identities" begin
    data = KC.LorentzianSpectralData(7 // 3, 4 // 5, 1 // 1)

    @test KC.spectral_energy(data) == 7 // 3
    @test KC.spectral_linewidth(data) == 4 // 5
    @test KC.spectral_residue(data) == 1 // 1
    @test KC.spectral_normalization(data) == 1 // 1
    @test KC.spectral_squared_weight(data) == 5 // 2
    @test KC.lorentzian_integrated_power(data, 3) == 75 // 8
    @test KC.lorentzian_spectral_value(data, 7 // 3) == 5 // 1

    @test_throws DomainError KC.LorentzianSpectralData(0 // 1, 0 // 1)
    @test_throws ArgumentError KC.lorentzian_integrated_power(data, 0)
end

@testset "retarded self-energy convention" begin
    data = @inferred KC.lorentzian_spectral_data_from_retarded_self_energy(
        3 // 2, Complex(1 // 4, -2 // 5); residue=5 // 6
    )
    @test KC.spectral_energy(data) == 7 // 4
    @test KC.spectral_linewidth(data) == 4 // 5
    @test KC.spectral_residue(data) == 5 // 6

    @test_throws DomainError KC.lorentzian_spectral_data_from_retarded_self_energy(
        1 // 1, Complex(0 // 1, 1 // 4)
    )
end

@testset "factorized repeated-line weight and width power counting" begin
    basis = KC.MomentumBasis(3)
    q1 = KC.basis_momentum(basis, 2)
    q2 = KC.basis_momentum(basis, 3)
    term = finite_width_term(
        basis,
        KC.KineticLine{Boson}[
            finite_width_line(finite_width_ϕ, q1),
            finite_width_line(finite_width_ϕ, q1),
            finite_width_line(finite_width_ϕ, q2),
        ],
        Val(3),
    )

    reduction = @inferred KC.lorentzian_spectral_reduction(term)
    @test reduction isa KC.LorentzianSpectralReduction{Boson}
    @test KC.spectral_reduction_kind(reduction) === KC.LorentzianResolved
    @test KC.spectral_reduction_resolved(reduction)
    @test KC.spectral_jacobian(reduction) == 1 // 1
    @test length(KC.spectral_factors(reduction)) == 2

    repeated = only(
        factor for
        factor in KC.spectral_factors(reduction) if KC.spectral_multiplicity(factor) == 2
    )
    ordinary = only(
        factor for
        factor in KC.spectral_factors(reduction) if KC.spectral_multiplicity(factor) == 1
    )
    @test KC.linewidth_exponent(repeated) == -1
    @test KC.linewidth_exponent(ordinary) == 0

    repeated_data = KC.LorentzianSpectralData(2 // 1, 4 // 5)
    ordinary_data = KC.LorentzianSpectralData(3 // 1, 3 // 2)
    model = KC.LorentzianSpectralModel(
        Dict(
            KC.spectral_line(repeated) => repeated_data,
            KC.spectral_line(ordinary) => ordinary_data,
        ),
    )
    @test @inferred(KC.evaluate_spectral_weight(reduction, model)) == 5 // 2

    counting = @inferred KC.width_aware_power_counting(
        KC.ParameterMonomial(:γ)^2, reduction
    )
    @test KC.parameters(counting) == KC.ParameterMonomial(:γ)^2
    powers = KC.linewidth_powers(counting)
    @test length(powers) == 1
    @test KC.spectral_line(only(powers)) == KC.spectral_line(repeated)
    @test KC.linewidth_exponent(only(powers)) == -1

    incomplete = KC.LorentzianSpectralModel(
        Dict(KC.spectral_line(ordinary) => ordinary_data)
    )
    @test_throws ArgumentError KC.evaluate_spectral_weight(reduction, incomplete)
end

@testset "one-residual-shell Lorentzian convolution" begin
    basis = KC.MomentumBasis(3)
    q1 = KC.basis_momentum(basis, 2)
    q2 = KC.basis_momentum(basis, 3)
    q12 = q1 + q2
    term = finite_width_term(
        basis,
        KC.KineticLine{Boson}[
            finite_width_line(finite_width_ϕ, q1),
            finite_width_line(finite_width_ϕ, q2),
            finite_width_line(finite_width_ϕ, q12),
        ],
        Val(3),
    )

    reduction = @inferred KC.lorentzian_convolution_reduction(term, finite_width_ϕ)
    @test reduction isa KC.LorentzianConvolutionReduction{Boson}
    @test KC.convolution_reduction_kind(reduction) === KC.LorentzianConvolutionResolved
    @test KC.convolution_reduction_resolved(reduction)
    @test KC.convolution_jacobian(reduction) == 1 // 1
    @test length(KC.convolution_pivot_lines(reduction)) == 2
    @test length(KC.convolution_coefficients(reduction)) == 2
    @test all(==(1 // 1), KC.convolution_coefficients(reduction))
    @test KC.convolution_external_coefficient(reduction) == 0 // 1

    lines = vcat(
        KC.convolution_pivot_lines(reduction), [KC.convolution_dependent_line(reduction)]
    )
    data = Dict(
        line => begin
            routed = KC.spectral_line_momentum(line)
            if routed == q1
                KC.LorentzianSpectralData(1 // 1, 1 // 1)
            elseif routed == q2
                KC.LorentzianSpectralData(2 // 1, 2 // 1)
            elseif routed == q12
                KC.LorentzianSpectralData(3 // 1, 3 // 1)
            else
                error("unexpected synthetic convolution momentum")
            end
        end for line in lines
    )
    model = KC.LorentzianSpectralModel(data)
    @test @inferred(KC.convolution_effective_linewidth(reduction, model)) == 6 // 1
    @test @inferred(KC.evaluate_spectral_convolution(reduction, model, 7 // 1)) == 2 // 3

    strict = @inferred KC.full_rank_frequency_reduction(term, finite_width_ϕ)
    shell, support_factor = KC.energy_shell(KC.convolution_energy_mismatch(reduction))
    support = KC.frequency_support(strict)
    @test support.shells == [shell]
    @test isempty(support.principal_values)
    @test KC.frequency_factor(strict) == KC.convolution_jacobian(reduction) * support_factor
end

@testset "finite-width analytic backends keep unsupported structures explicit" begin
    basis = KC.MomentumBasis(3)
    q1 = KC.basis_momentum(basis, 2)
    q2 = KC.basis_momentum(basis, 3)
    q12 = q1 + q2

    convolution = finite_width_term(
        basis,
        KC.KineticLine{Boson}[
            finite_width_line(finite_width_ϕ, q1),
            finite_width_line(finite_width_ϕ, q2),
            finite_width_line(finite_width_ϕ, q12),
        ],
        Val(3),
    )
    reduction = @inferred KC.lorentzian_spectral_reduction(convolution)
    @test reduction isa KC.LorentzianSpectralReduction{Boson}
    @test KC.spectral_reduction_kind(reduction) === KC.LorentzianUnsupported
    @test !KC.spectral_reduction_resolved(reduction)
    @test_throws ArgumentError KC.spectral_jacobian(reduction)
    @test_throws ArgumentError KC.spectral_factors(reduction)

    repeated = finite_width_term(
        basis,
        KC.KineticLine{Boson}[
            finite_width_line(finite_width_ϕ, q1),
            finite_width_line(finite_width_ϕ, q1),
            finite_width_line(finite_width_ϕ, q2),
        ],
        Val(3),
    )
    convolution_reduction = @inferred KC.lorentzian_convolution_reduction(
        repeated, finite_width_ϕ
    )
    @test convolution_reduction isa KC.LorentzianConvolutionReduction{Boson}
    @test KC.convolution_reduction_kind(convolution_reduction) ===
        KC.LorentzianConvolutionUnsupported
    @test !KC.convolution_reduction_resolved(convolution_reduction)
end
