using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields finite_width_variation_ϕ::Boson

function variation_line(
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

function variation_term(
    basis::KC.MomentumBasis, lines::Vector{KC.KineticLine{S}}, ::Val{E}
) where {S<:KC.Statistics,E}
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

@testset "integrated Lorentzian power variation" begin
    data = KC.LorentzianSpectralData(7 // 3, 4 // 5, 5 // 6)
    variation = KC.LorentzianSpectralDataVariation(11 // 7, 2 // 3, 3 // 5)

    @test KC.spectral_energy_variation(variation) == 11 // 7
    @test KC.spectral_linewidth_variation(variation) == 2 // 3
    @test KC.spectral_residue_variation(variation) == 3 // 5
    @test @inferred(KC.lorentzian_integrated_power_variation(data, variation, 1)) == 3 // 5

    Z = KC.spectral_residue(data)
    Γ = KC.spectral_linewidth(data)
    δZ = KC.spectral_residue_variation(variation)
    δΓ = KC.spectral_linewidth_variation(variation)
    expected = 4 * Z * δZ / Γ - 2 * Z^2 * δΓ / Γ^2
    @test @inferred(KC.lorentzian_integrated_power_variation(data, variation, 2)) ==
        expected

    energy_only = KC.LorentzianSpectralDataVariation(9 // 4, 0 // 1, 0 // 1)
    @test iszero(KC.lorentzian_integrated_power_variation(data, energy_only, 3))
end

@testset "factorized repeated-line spectral variation" begin
    basis = KC.MomentumBasis(3)
    q1 = KC.basis_momentum(basis, 2)
    q2 = KC.basis_momentum(basis, 3)
    term = variation_term(
        basis,
        KC.KineticLine{Boson}[
            variation_line(finite_width_variation_ϕ, q1),
            variation_line(finite_width_variation_ϕ, q1),
            variation_line(finite_width_variation_ϕ, q2),
        ],
        Val(3),
    )
    reduction = KC.lorentzian_spectral_reduction(term)
    @test KC.spectral_reduction_resolved(reduction)

    factors = KC.spectral_factors(reduction)
    repeated = only(factor for factor in factors if KC.spectral_multiplicity(factor) == 2)
    ordinary = only(factor for factor in factors if KC.spectral_multiplicity(factor) == 1)
    repeated_line = KC.spectral_line(repeated)
    ordinary_line = KC.spectral_line(ordinary)

    repeated_data = KC.LorentzianSpectralData(2 // 1, 4 // 5, 5 // 6)
    ordinary_data = KC.LorentzianSpectralData(3 // 1, 3 // 2, 7 // 8)
    repeated_variation = KC.LorentzianSpectralDataVariation(1 // 9, 2 // 3, 3 // 5)
    ordinary_variation = KC.LorentzianSpectralDataVariation(-2 // 7, 5 // 11, -1 // 4)
    model = KC.LorentzianSpectralModel(
        Dict(repeated_line => repeated_data, ordinary_line => ordinary_data)
    )
    model_variation = KC.LorentzianSpectralModelVariation(
        Dict(repeated_line => repeated_variation, ordinary_line => ordinary_variation)
    )

    B2 = KC.lorentzian_integrated_power(repeated_data, 2)
    δB2 = KC.lorentzian_integrated_power_variation(repeated_data, repeated_variation, 2)
    Z2 = KC.spectral_residue(ordinary_data)
    δZ2 = KC.spectral_residue_variation(ordinary_variation)
    expected = KC.spectral_jacobian(reduction) * (δB2 * Z2 + B2 * δZ2)
    @test @inferred(
        KC.evaluate_spectral_weight_variation(reduction, model, model_variation)
    ) == expected

    missing = KC.LorentzianSpectralModelVariation(Dict(repeated_line => repeated_variation))
    @test_throws ArgumentError KC.evaluate_spectral_weight_variation(
        reduction, model, missing
    )
end

@testset "residual convolution spectral variation" begin
    basis = KC.MomentumBasis(3)
    q1 = KC.basis_momentum(basis, 2)
    q2 = KC.basis_momentum(basis, 3)
    q12 = q1 + q2
    term = variation_term(
        basis,
        KC.KineticLine{Boson}[
            variation_line(finite_width_variation_ϕ, q1),
            variation_line(finite_width_variation_ϕ, q2),
            variation_line(finite_width_variation_ϕ, q12),
        ],
        Val(3),
    )
    reduction = KC.lorentzian_convolution_reduction(term, finite_width_variation_ϕ)
    @test KC.convolution_reduction_resolved(reduction)
    @test all(==(1 // 1), KC.convolution_coefficients(reduction))
    @test KC.convolution_external_coefficient(reduction) == 0 // 1

    lines = vcat(
        KC.convolution_pivot_lines(reduction), [KC.convolution_dependent_line(reduction)]
    )
    data = Dict(
        line => begin
            routed = KC.spectral_line_momentum(line)
            if routed == q1
                KC.LorentzianSpectralData(1 // 1, 1 // 1, 2 // 1)
            elseif routed == q2
                KC.LorentzianSpectralData(2 // 1, 2 // 1, 3 // 1)
            elseif routed == q12
                KC.LorentzianSpectralData(4 // 1, 3 // 1, 5 // 1)
            else
                error("unexpected synthetic momentum")
            end
        end for line in lines
    )
    variations = Dict(
        line => begin
            routed = KC.spectral_line_momentum(line)
            if routed == q1
                KC.LorentzianSpectralDataVariation(1 // 7, 1 // 5, 1 // 3)
            elseif routed == q2
                KC.LorentzianSpectralDataVariation(-2 // 7, -1 // 6, 2 // 5)
            elseif routed == q12
                KC.LorentzianSpectralDataVariation(3 // 8, 1 // 4, -1 // 2)
            else
                error("unexpected synthetic momentum")
            end
        end for line in lines
    )
    model = KC.LorentzianSpectralModel(data)
    model_variation = KC.LorentzianSpectralModelVariation(variations)

    @test KC.evaluate_spectral_convolution(reduction, model, 0 // 1) == 18 // 1
    @test @inferred(
        KC.evaluate_spectral_convolution_variation(
            reduction, model, model_variation, 0 // 1
        )
    ) == 739 // 700

    zero_variation = KC.LorentzianSpectralModelVariation(
        Dict(line => KC.LorentzianSpectralDataVariation(0, 0, 0) for line in lines)
    )
    @test iszero(
        KC.evaluate_spectral_convolution_variation(reduction, model, zero_variation, 0 // 1)
    )
end

@testset "spectral variation API is qualified-public" begin
    for name in (
        :LorentzianSpectralDataVariation,
        :LorentzianSpectralModelVariation,
        :spectral_data_variation,
        :spectral_energy_variation,
        :spectral_linewidth_variation,
        :spectral_residue_variation,
        :lorentzian_integrated_power_variation,
        :evaluate_spectral_weight_variation,
        :evaluate_spectral_convolution_variation,
    )
        if VERSION >= v"1.11"
            @test !Base.isexported(KC, name)
            @test Base.ispublic(KC, name)
        else
            @test !(name in Base.names(KC; all=false, imported=false))
        end
    end
end
