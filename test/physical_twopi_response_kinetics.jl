using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields response_kinetics_ψ::Boson response_kinetics_χ::Boson

function response_kinetics_hs_interaction()
    ψc = response_kinetics_ψ[Classical]
    ψq = response_kinetics_ψ[Quantum]
    χc = response_kinetics_χ[Classical]
    χq = response_kinetics_χ[Quantum]

    forward = (1 // 2) * ψc^2 * bar(χq) + ψc * ψq * bar(χc) + (1 // 2) * ψq^2 * bar(χq)
    return ChargedInteractionLagrangian(
        -im * (forward + bar(forward)),
        response_kinetics_ψ => 1,
        response_kinetics_χ => 2;
        parameter=:h,
    )
end

function response_kinetics_components(object)
    return (
        KC.keldysh_component(object),
        KC.retarded_component(object),
        KC.advanced_component(object),
    )
end

@testset "KC-native exact regular HS response" begin
    g = 3 // 2
    γ = 2 // 5
    ΩK = complex(5 // 13, -2 // 7)
    ΩR = complex(2 // 7, 5 // 11)
    ΩA = complex(-3 // 10, 4 // 9)

    response = KC._two_body_loss_regular_response_components(g, γ, ΩK, ΩR, ΩA)
    λ = complex(g, -γ)
    λbar = conj(λ)

    DR = -im * λ + response.retarded
    DA = -im * λbar + response.advanced
    DK = response.full_keldysh

    # Exact component Dyson identities for D = D0 + D0 Ω D.
    @test DR == -im * λ - im * λ * ΩR * DR
    @test DA == -im * λbar - im * λbar * ΩA * DA
    @test DK == -2 * γ - im * λ * ΩR * DK - 2 * γ * ΩA * DA - im * λ * ΩK * DA
    @test response.keldysh == DK + 2 * γ

    # The linear-in-Ω coefficients are the frozen C3 native response table. Since Dreg has no
    # constant term, denominator corrections start only at quadratic order in Ω.
    @test -(λ^2) == -g^2 + 2im * g * γ + γ^2
    @test -(λbar^2) == -g^2 - 2im * g * γ + γ^2
    @test -(λ * λbar) == -g^2 - γ^2
    @test 2im * γ * λ == 2im * g * γ + 2 * γ^2
    @test 2im * γ * λbar == 2im * g * γ - 2 * γ^2
end

@testset "response-aware homogeneous 2PI kinetic lowering" begin
    Γ2 = TwoPIEffectiveAction(response_kinetics_hs_interaction(), Val(2), Val(3))
    carrier = KC.twopi_composite_fourier_self_energy(
        Γ2, response_kinetics_ψ, response_kinetics_χ
    )

    # The ordinary path remains forbidden: χ must never become an ordinary Wigner/kinetic line.
    @test_throws ArgumentError wigner_transform(carrier)
    @test_throws ArgumentError KC.response_wigner_transform(carrier; gradient_order=Val(1))

    ΣW = KC.response_wigner_transform(carrier)
    @test KC.target_family(ΣW) === response_kinetics_ψ
    @test KC.statistics(ΣW) === Boson
    @test KC.gradient_order(ΣW) == Val(0)

    ΩW = KC.response_polarization(ΣW)
    @test KC.response_physical_family(ΩW) === response_kinetics_ψ
    @test KC.response_coherent_parameter(ΩW) === :g
    @test KC.response_loss_parameter(ΩW) === :γ

    # The nested polarization has its own local Wigner basis. Its external pair momentum is the
    # first local basis variable and its internal lines are all physical ψ lines.
    for component in response_kinetics_components(ΩW)
        @test !isempty(component)
        for (graph, contributions) in component
            @test KC.external_wigner_momentum(graph) == KC.MomentumVariable(1)
            @test all(
                edge -> KC.field_family(edge.out) === response_kinetics_ψ,
                KC.contractions(KC.coordinate_diagram(graph)),
            )
            @test !isempty(contributions)
        end
    end

    # The outer Wigner carrier has eliminated χ completely. Only the response component and its
    # routed pair momentum survive as the binding from the outer loop basis into the local pair
    # polarization basis.
    for component in response_kinetics_components(ΣW)
        @test !isempty(component)
        for (graph, contributions) in component
            @test length(KC.physical_contractions(graph)) == 1
            @test all(
                edge -> KC.field_family(edge.out) === response_kinetics_ψ,
                KC.physical_contractions(graph),
            )
            @test KC.response_component(graph) in (
                KC.PropagatorType.Keldysh,
                KC.PropagatorType.Retarded,
                KC.PropagatorType.Advanced,
            )
            @test length(KC.response_momentum(graph)) == length(KC.momentum_basis(graph))
            @test !isempty(contributions)
        end
    end

    ΣK = kinetic_expression(ΣW)
    @test KC.target_family(ΣK) === response_kinetics_ψ
    @test KC.statistics(ΣK) === Boson
    @test KC.gradient_order(ΣK) == Val(0)

    ΩKinetic = KC.response_polarization(ΣK)
    @test KC.response_physical_family(ΩKinetic) === response_kinetics_ψ
    for component in response_kinetics_components(ΩKinetic)
        @test !isempty(component)
        for (term, _) in component
            @test KC.external_wigner_momentum(term) == KC.MomentumVariable(1)
            @test all(line -> line.family === response_kinetics_ψ, KC.kinetic_lines(term))
        end
    end

    for component in response_kinetics_components(ΣK)
        @test !isempty(component)
        for (term, _) in component
            @test length(KC.kinetic_lines(term)) == 1
            @test all(line -> line.family === response_kinetics_ψ, KC.kinetic_lines(term))
            @test KC.response_component(term) in (
                KC.PropagatorType.Keldysh,
                KC.PropagatorType.Retarded,
                KC.PropagatorType.Advanced,
            )
            @test length(KC.response_momentum(term)) == length(KC.momentum_basis(term))
        end
    end

    collision = off_shell_collision_expression(ΣK)
    @test KC.target_family(collision) === response_kinetics_ψ
    @test KC.statistics(collision) === Boson
    @test KC.gradient_order(collision) == Val(0)
    @test !isempty(KC.collision_offset(collision))
    @test !isempty(KC.collision_distribution_coefficient(collision))
    @test KC.response_physical_family(KC.response_polarization(collision)) ===
        response_kinetics_ψ

    # No quasiparticle shell is assumed for the exact resummed pair continuum.
    @test_throws ArgumentError spectral_dispersive_collision(collision)
end
