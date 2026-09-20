using KeldyshContraction, Test
import KeldyshContraction as KC

@testset "exact KC-native regular HS response" begin
    g = 3 // 2
    γ = 1 // 5
    ΩK = 2 // 7
    ΩR = 1 // 11
    ΩA = 3 // 13

    response = @inferred KC._two_body_loss_regular_response(g, γ, ΩK, ΩR, ΩA)

    λ = g - im * γ
    λbar = g + im * γ
    D0 = Complex{Rational{Int}}[
        -2γ (-im * λ)
        (-im * λbar) 0
    ]
    Ω = Complex{Rational{Int}}[0 ΩA; ΩR ΩK]
    Dreg = Complex{Rational{Int}}[
        response.keldysh response.retarded
        response.advanced 0
    ]
    D = D0 + Dreg

    @test D == D0 + D0 * Ω * D

    linear = KC._two_body_loss_regular_response_linear(g, γ, ΩK, ΩR, ΩA)
    sectors = KC._two_body_loss_regular_response_linear_sectors(ΩK, ΩR, ΩA)
    reconstructed = (
        keldysh=g^2 * sectors.g2.keldysh +
                g * γ * sectors.gγ.keldysh +
                γ^2 * sectors.γ2.keldysh,
        retarded=g^2 * sectors.g2.retarded +
                 g * γ * sectors.gγ.retarded +
                 γ^2 * sectors.γ2.retarded,
        advanced=g^2 * sectors.g2.advanced +
                 g * γ * sectors.gγ.advanced +
                 γ^2 * sectors.γ2.advanced,
    )
    @test linear == reconstructed

    @test sectors.g2 == (keldysh=(-ΩK), retarded=(-ΩR), advanced=(-ΩA))
    @test sectors.gγ == (keldysh=2im * (ΩR + ΩA), retarded=2im * ΩR, advanced=-2im * ΩA)
    @test sectors.γ2 == (keldysh=(-ΩK + 2ΩR - 2ΩA), retarded=ΩR, advanced=ΩA)
end

@qfields response_wigner_ψ::Boson response_wigner_χ::Boson

function response_wigner_hs_interaction()
    ψc = response_wigner_ψ[Classical]
    ψq = response_wigner_ψ[Quantum]
    χc = response_wigner_χ[Classical]
    χq = response_wigner_χ[Quantum]

    forward = (1 // 2) * ψc^2 * bar(χq) + ψc * ψq * bar(χc) + (1 // 2) * ψq^2 * bar(χq)
    return ChargedInteractionLagrangian(
        -im * (forward + bar(forward)),
        response_wigner_ψ => 1,
        response_wigner_χ => 2;
        parameter=:h,
    )
end

@testset "response-aware homogeneous Wigner carrier" begin
    Γ2 = @inferred TwoPIEffectiveAction(response_wigner_hs_interaction(), Val(2), Val(3))
    fourier = @inferred KC.twopi_composite_fourier_self_energy(
        Γ2, response_wigner_ψ, response_wigner_χ; coherent_parameter=:g, loss_parameter=:γ
    )
    wigner = @inferred KC.response_aware_wigner_transform(fourier)

    @test KC.target_family(wigner) === response_wigner_ψ
    @test KC.response_family(wigner) === response_wigner_χ
    @test KC.response_polarization(wigner) === KC.response_polarization(fourier)
    @test KC.gradient_order(wigner) == Val(0)
    @test KC.wigner_context(wigner) isa KC.HomogeneousWignerContext

    for component in (
        KC.keldysh_component(wigner),
        KC.retarded_component(wigner),
        KC.advanced_component(wigner),
    )
        @test !isempty(component)
        for (graph, contributions) in component
            @test KC.response_family(graph) === response_wigner_χ
            @test KC.response_component(graph) in (
                KC.PropagatorType.Keldysh,
                KC.PropagatorType.Retarded,
                KC.PropagatorType.Advanced,
            )
            @test KC.response_momentum(graph) == KC.response_momentum(graph.composite)
            @test KC.external_wigner_momentum(graph) == KC.momentum_basis(graph)[1]
            @test all(
                edge -> KC.field_family(edge.out) === response_wigner_ψ,
                KC.physical_contractions(graph),
            )
            @test !isempty(contributions)
        end
    end

    @test_throws ArgumentError KC.response_aware_wigner_transform(
        fourier; gradient_order=Val(1)
    )
end
