using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields composite_response_ψ::Boson composite_response_χ::Boson

function composite_response_hs_interaction()
    ψc = composite_response_ψ[Classical]
    ψq = composite_response_ψ[Quantum]
    χc = composite_response_χ[Classical]
    χq = composite_response_χ[Quantum]

    # KC-native symmetric Keldysh rotation certified in C3.
    forward =
        (1 // 2) * ψc^2 * bar(χq) + ψc * ψq * bar(χc) + (1 // 2) * ψq^2 * bar(χq)
    return ChargedInteractionLagrangian(
        -im * (forward + bar(forward)),
        composite_response_ψ => 1,
        composite_response_χ => 2;
        parameter=:h,
    )
end

@testset "anchored composite HS response carrier" begin
    Γ2 = @inferred TwoPIEffectiveAction(composite_response_hs_interaction(), Val(2), Val(3))
    Σ = @inferred KC.twopi_composite_fourier_self_energy(
        Γ2,
        composite_response_ψ,
        composite_response_χ;
        coherent_parameter=:g,
        loss_parameter=:γ,
    )

    @test KC.target_family(Σ) === composite_response_ψ
    @test KC.statistics(Σ) === Boson
    @test KC.response_family(Σ) === composite_response_χ
    @test KC.response_coherent_parameter(Σ) === :g
    @test KC.response_loss_parameter(Σ) === :γ

    Ω = KC.response_polarization(Σ)
    @test KC.target_family(Ω) === composite_response_χ
    for component in (
        KC.keldysh_component(Ω), KC.retarded_component(Ω), KC.advanced_component(Ω)
    )
        @test !isempty(component)
        for (graph, contributions) in component
            @test graph.external_count == 1
            @test graph.loop_count == 1
            @test length(KC.contractions(graph.coordinate)) == 2
            @test all(
                edge -> KC.field_family(edge.out) === composite_response_ψ,
                KC.contractions(graph.coordinate),
            )
            @test !isempty(contributions)
        end
    end

    components = (
        KC.keldysh_component(Σ), KC.retarded_component(Σ), KC.advanced_component(Σ)
    )
    for component in components
        @test !isempty(component)
        for (graph, contributions) in component
            @test KC.external_momentum_count(graph) == 1
            @test KC.loop_momentum_count(graph) == 1
            @test length(KC.contractions(KC.coordinate_diagram(graph))) == 2
            @test KC.response_family(graph) === composite_response_χ
            @test KC.response_component(graph) in (
                KC.PropagatorType.Keldysh,
                KC.PropagatorType.Retarded,
                KC.PropagatorType.Advanced,
            )
            @test KC.response_momentum(graph) ==
                KC.edge_momenta(graph)[KC.response_edge_index(graph)]

            physical = KC.physical_contractions(graph)
            @test length(physical) == 1
            @test all(
                edge -> KC.field_family(edge.out) === composite_response_ψ, physical
            )
            @test !isempty(contributions)
        end
    end

    error = try
        wigner_transform(Σ)
        nothing
    catch exception
        exception
    end
    @test error isa ArgumentError
    @test occursin("eliminated into physical ψ content", sprint(showerror, error))
end

@testset "composite response constructor rejects invalid provenance" begin
    Γ2 = TwoPIEffectiveAction(composite_response_hs_interaction(), Val(2), Val(3))

    @test_throws ArgumentError KC.twopi_composite_fourier_self_energy(
        Γ2, composite_response_ψ, composite_response_ψ
    )

    polarization = KC.twopi_fourier_self_energy(Γ2, composite_response_χ)
    @test_throws ArgumentError KC.TwoBodyLossHSResponse(
        polarization,
        composite_response_χ;
        coherent_parameter=:g,
        loss_parameter=:g,
    )
end
