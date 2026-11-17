using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields response_collision_ψ::Boson response_collision_χ::Boson

function response_collision_hs_interaction()
    ψc = response_collision_ψ[Classical]
    ψq = response_collision_ψ[Quantum]
    χc = response_collision_χ[Classical]
    χq = response_collision_χ[Quantum]

    forward = (1 // 2) * ψc^2 * bar(χq) + ψc * ψq * bar(χc) + (1 // 2) * ψq^2 * bar(χq)
    return ChargedInteractionLagrangian(
        -im * (forward + bar(forward)),
        response_collision_ψ => 1,
        response_collision_χ => 2;
        parameter=:h,
    )
end

@testset "response-aware physical kinetic lowering" begin
    Γ2 = @inferred TwoPIEffectiveAction(response_collision_hs_interaction(), Val(2), Val(3))
    composite = @inferred KC.twopi_composite_fourier_self_energy(
        Γ2,
        response_collision_ψ,
        response_collision_χ;
        coherent_parameter=:g,
        loss_parameter=:γ,
    )
    wigner = @inferred KC.response_aware_wigner_transform(composite)
    kinetic = @inferred kinetic_expression(wigner)

    @test KC.target_family(kinetic) === response_collision_ψ
    @test KC.response_family(kinetic) === response_collision_χ
    @test KC.gradient_order(kinetic) == Val(0)

    # The nested χ-target polarization is allowed to retain χ as target provenance, but every
    # actual kinetic propagator line inside it is physical ψ.
    polarization = KC.response_polarization(kinetic)
    @test KC.target_family(polarization) === response_collision_χ
    for expression in (polarization.keldysh, polarization.retarded, polarization.advanced)
        for (term, _) in expression
            @test all(line -> line.family === response_collision_ψ, KC.kinetic_lines(term))
        end
    end

    for expression in (
        KC.keldysh_component(kinetic),
        KC.retarded_component(kinetic),
        KC.advanced_component(kinetic),
    )
        @test !isempty(expression)
        for (term, _) in expression
            @test length(KC.kinetic_lines(term)) == 1
            @test all(line -> line.family === response_collision_ψ, KC.kinetic_lines(term))
            @test KC.response_component(term) in (
                KC.PropagatorType.Keldysh,
                KC.PropagatorType.Retarded,
                KC.PropagatorType.Advanced,
            )
        end
    end
end

@testset "response-aware off-shell KB collision" begin
    Γ2 = TwoPIEffectiveAction(response_collision_hs_interaction(), Val(2), Val(3))
    composite = KC.twopi_composite_fourier_self_energy(
        Γ2,
        response_collision_ψ,
        response_collision_χ;
        coherent_parameter=:g,
        loss_parameter=:γ,
    )
    kinetic = kinetic_expression(KC.response_aware_wigner_transform(composite))
    collision = @inferred off_shell_collision_expression(kinetic)

    @test KC.target_family(collision) === response_collision_ψ
    @test KC.response_family(collision) === response_collision_χ
    @test KC.response_polarization(collision) === KC.response_polarization(kinetic)
    @test KC.gradient_order(collision) == Val(0)

    offset = KC.collision_offset(collision)
    distribution = KC.collision_distribution_coefficient(collision)
    @test offset.terms == (im * kinetic.keldysh).terms
    @test distribution.terms == (-im * (kinetic.retarded - kinetic.advanced)).terms

    # No χ line can enter either affine collision coefficient.
    for expression in (offset, distribution)
        for (term, _) in expression
            @test all(line -> line.family === response_collision_ψ, KC.kinetic_lines(term))
        end
    end
end
