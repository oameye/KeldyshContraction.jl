using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields external_projection_ϕ::Boson

function external_projection_line(
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

function external_projection_term(
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

function generated_external_projection_collision()
    c = external_projection_ϕ[Classical]
    q = external_projection_ϕ[Quantum]
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus
    loss =
        (1 // 2) *
        im *
        (
            bar(c) * bar(q) * (c(minus)^2 + q(minus)^2) -
            c(plus) * q(plus) * (bar(c)^2 + bar(q)^2) +
            2 * bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
        )
    L = InteractionLagrangian(loss, :γ)
    G = DressedPropagator(L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false)
    ΣF = SelfEnergy(fourier_transform(G))
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    kinetic = kinetic_expression(ΣW)
    return spectral_dispersive_collision(off_shell_collision_expression(kinetic))
end

function generated_external_projection_regular_sector(collision)
    target = KC.target_family(collision)
    for expression in
        (KC.collision_offset(collision), KC.collision_distribution_coefficient(collision))
        for (term, coefficient) in expression
            iszero(coefficient) && continue
            kinds = KC.spectral_dispersive_kinds(term)
            count(==(KC.CollisionSpectral), kinds) == 3 || continue
            lines = KC.kinetic_lines(term.carrier)
            identities = [
                KC.SpectralLineIdentity{Boson}(line.family, KC.momentum(line)) for
                (line, kind) in zip(lines, kinds) if kind === KC.CollisionSpectral
            ]
            length(unique(identities)) == length(identities) || continue
            reduction = KC.lorentzian_convolution_reduction(term, target)
            KC.convolution_reduction_resolved(reduction) || continue
            return term, reduction
        end
    end
    return error("generated loss collision contains no resolved regular convolution sector")
end

function external_projection_external_line(term, target)
    basis = KC.momentum_basis(term)
    external_index = KC._external_frequency_basis_index(term)
    return KC.SpectralLineIdentity{Boson}(target, KC.basis_momentum(basis, external_index))
end

@testset "external Lorentzian projection adds the full spectral linewidth" begin
    basis = KC.MomentumBasis(3)
    k = KC.basis_momentum(basis, 1)
    q1 = KC.basis_momentum(basis, 2)
    q2 = KC.basis_momentum(basis, 3)
    dependent = k - q1 - q2
    term = external_projection_term(
        basis,
        KC.KineticLine{Boson}[
            external_projection_line(external_projection_ϕ, q1),
            external_projection_line(external_projection_ϕ, q2),
            external_projection_line(external_projection_ϕ, dependent),
        ],
        Val(3),
    )
    reduction = @inferred KC.lorentzian_convolution_reduction(term, external_projection_ϕ)
    @test KC.convolution_reduction_resolved(reduction)
    @test abs(KC.convolution_external_coefficient(reduction)) == 1 // 1
    @test all(
        coefficient -> abs(coefficient) == 1 // 1, KC.convolution_coefficients(reduction)
    )
    @test KC.convolution_jacobian(reduction) == 1 // 1

    external = KC.SpectralLineIdentity{Boson}(external_projection_ϕ, k)
    internal = vcat(
        KC.convolution_pivot_lines(reduction), [KC.convolution_dependent_line(reduction)]
    )
    widths = Dict(internal[1] => 1 // 1, internal[2] => 3 // 1, internal[3] => 4 // 1)
    data = Dict(
        line => KC.LorentzianSpectralData(0 // 1, width, 1 // 1) for (line, width) in widths
    )
    data[external] = KC.LorentzianSpectralData(0 // 1, 2 // 1, 1 // 1)
    model = KC.LorentzianSpectralModel(data)

    @test KC.convolution_effective_linewidth(reduction, model) == 8 // 1
    @test @inferred(
        KC.external_spectral_projection_linewidth(reduction, model, external)
    ) == 10 // 1
    @test @inferred(KC.external_spectral_projection_mismatch(reduction, model, external)) ==
        0 // 1
    @test @inferred(KC.evaluate_spectral_convolution(reduction, model, 0 // 1)) == 1 // 2
    @test @inferred(KC.evaluate_external_spectral_projection(reduction, model, external)) ==
        2 // 5

    variation = KC.LorentzianSpectralModelVariation(
        Dict(
            line => KC.LorentzianSpectralDataVariation(0 // 1, 1 // 1, 0 // 1) for
            line in keys(data)
        ),
    )
    @test @inferred(
        KC.evaluate_external_spectral_projection_variation(
            reduction, model, variation, external
        )
    ) == -4 // 25

    incomplete = KC.LorentzianSpectralModel(Dict(line => data[line] for line in internal))
    @test_throws ArgumentError KC.evaluate_external_spectral_projection(
        reduction, incomplete, external
    )
end

@testset "generated γ² sunset receives the fourth spectral linewidth" begin
    collision = generated_external_projection_collision()
    term, reduction = generated_external_projection_regular_sector(collision)
    target = KC.target_family(collision)
    external = external_projection_external_line(term, target)
    β = KC.convolution_external_coefficient(reduction)

    @test abs(β) == 1 // 1
    @test all(
        coefficient -> abs(coefficient) == 1 // 1, KC.convolution_coefficients(reduction)
    )

    internal = vcat(
        KC.convolution_pivot_lines(reduction), [KC.convolution_dependent_line(reduction)]
    )
    lines = unique(vcat(internal, [external]))
    @test length(lines) == 4
    model = KC.LorentzianSpectralModel(
        Dict(line => KC.LorentzianSpectralData(0 // 1, 1 // 1, 1 // 1) for line in lines)
    )
    variation = KC.LorentzianSpectralModelVariation(
        Dict(
            line => KC.LorentzianSpectralDataVariation(0 // 1, 1 // 1, 0 // 1) for
            line in lines
        ),
    )

    internal_width = KC.convolution_effective_linewidth(reduction, model)
    full_width = KC.external_spectral_projection_linewidth(reduction, model, external)
    @test internal_width == 3 // 1
    @test full_width == 4 // 1
    @test KC.external_spectral_projection_mismatch(reduction, model, external) == 0 // 1

    jacobian = KC.convolution_jacobian(reduction)
    @test KC.evaluate_external_spectral_projection(reduction, model, external) == jacobian
    @test KC.evaluate_external_spectral_projection_variation(
        reduction, model, variation, external
    ) == -jacobian
end
