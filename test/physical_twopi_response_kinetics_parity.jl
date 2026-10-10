using KeldyshContraction, Test
import KeldyshContraction as KC

# Reuse the frozen C4 anchored weak-expansion oracle. The file also reruns its own C4 testsets in
# this isolated test sandbox; that is intentional here because C5 must certify against the exact
# same endpoint-preserving expansion rather than introduce a second attachment convention.
include(joinpath(@__DIR__, "physical_twopi_composite_response_parity.jl"))

function response_kinetics_component(object, kind::Symbol)
    kind === :K && return KC.keldysh_component(object)
    kind === :R && return KC.retarded_component(object)
    kind === :A && return KC.advanced_component(object)
    return error("unexpected response-aware component")
end

function response_kinetics_assert_c5_carrier(c5, carrier)
    expected = kinetic_expression(KC.response_wigner_transform(carrier))

    @test KC.target_family(c5) === composite_parity_ψ
    @test KC.target_family(expected) === composite_parity_ψ
    @test KC.statistics(c5) === Boson
    @test KC.gradient_order(c5) == Val(0)

    for kind in (:K, :R, :A)
        actual_component = response_kinetics_component(c5, kind)
        expected_component = response_kinetics_component(expected, kind)
        @test actual_component.terms == expected_component.terms
        for (term, _) in actual_component
            @test all(
                line -> line.family === composite_parity_ψ,
                KC.kinetic_lines(term),
            )
            @test KC.response_component(term) in (
                KC.PropagatorType.Keldysh,
                KC.PropagatorType.Retarded,
                KC.PropagatorType.Advanced,
            )
        end
    end

    actual_response = KC.response_polarization(c5)
    expected_response = KC.response_polarization(expected)
    @test KC.response_physical_family(actual_response) === composite_parity_ψ
    @test KC.response_physical_family(expected_response) === composite_parity_ψ
    for kind in (:K, :R, :A)
        @test response_kinetics_component(actual_response, kind) ==
            response_kinetics_component(expected_response, kind)
    end
    return nothing
end

"""
Test-only controlled weak expansion of the C5 response-aware object.

C5 has already eliminated the auxiliary χ edge before statistical lowering, so it deliberately no
longer carries the endpoint identity needed to attach the pair bubble as a finite graph. The weak
oracle therefore returns to the original Γ₂ cut provenance certified in C3/C4, but only after
checking that the supplied C5 object is exactly the response-aware lowering of that same carrier.
"""
function response_kinetics_weak_expansion(Γ2, c5, sector::Symbol, parameter)
    carrier = KC.twopi_composite_fourier_self_energy(
        Γ2, composite_parity_ψ, composite_parity_χ
    )
    response_kinetics_assert_c5_carrier(c5, carrier)
    return composite_parity_expanded_self_energy(Γ2, carrier, sector, parameter)
end

@testset "C5 controlled weak response returns to frozen C3 collision content" begin
    Γ2 = TwoPIEffectiveAction(composite_parity_hs_interaction(), Val(2), Val(3))
    carrier = KC.twopi_composite_fourier_self_energy(
        Γ2, composite_parity_ψ, composite_parity_χ
    )
    c5 = kinetic_expression(KC.response_wigner_transform(carrier))

    sectors = (
        (:g2, KC.ParameterMonomial(:g)^2),
        (:gγ, KC.ParameterMonomial(:g) * KC.ParameterMonomial(:γ)),
        (:γ2, KC.ParameterMonomial(:γ)^2),
    )

    expanded = Dict{Symbol,Any}()
    direct = Dict{Symbol,Any}()
    for (sector, parameter) in sectors
        expanded[sector] = response_kinetics_weak_expansion(Γ2, c5, sector, parameter)
        direct[sector] = composite_parity_microscopic_skeleton(parameter)
    end

    # The elastic sector is already identical at the raw Fourier boundary.
    @test expanded[:g2].keldysh == direct[:g2].keldysh
    @test expanded[:g2].retarded == direct[:g2].retarded
    @test expanded[:g2].advanced == direct[:g2].advanced

    # Mixed coherent-loss content agrees after the same physical frequency canonicalization used
    # by frozen C3. No auxiliary statistical species survives the C5 round trip.
    mixed_expanded = composite_parity_canonical_collision(expanded[:gγ])
    mixed_direct = composite_parity_canonical_collision(direct[:gγ])
    @test isempty(KC.shifted_frequency_terms(mixed_expanded))
    @test isempty(KC.shifted_frequency_terms(mixed_direct))
    @test KC.canonical_frequency_expressions(mixed_expanded) ==
        KC.canonical_frequency_expressions(mixed_direct)
    @test all(==(composite_parity_ψ), composite_parity_statistical_families(mixed_expanded))

    # The γ² 2PI skeleton has no finite-Trotter shifted sector. Its causal boundary reduction must
    # remain exactly the frozen C3 result after passing through the C5 response-aware object.
    gamma2_expanded = composite_parity_canonical_collision(expanded[:γ2])
    gamma2_direct = composite_parity_canonical_collision(direct[:γ2])
    @test isempty(KC.shifted_frequency_terms(gamma2_expanded))
    @test isempty(KC.shifted_frequency_terms(gamma2_direct))
    @test all(
        ==(composite_parity_ψ), composite_parity_statistical_families(gamma2_expanded)
    )

    expanded_reduction = KC.reduce_canonical_trotter_frequencies(gamma2_expanded)
    direct_reduction = KC.reduce_canonical_trotter_frequencies(gamma2_direct)
    @test isempty(KC.unresolved_trotter_states(expanded_reduction))
    @test isempty(KC.unresolved_trotter_states(direct_reduction))
    @test KC.canonical_trotter_constants(expanded_reduction) ==
        KC.canonical_trotter_constants(direct_reduction)
    @test KC.canonical_trotter_boundary_expressions(expanded_reduction) ==
        KC.canonical_trotter_boundary_expressions(direct_reduction)
end
