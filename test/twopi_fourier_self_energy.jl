using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields twopi_fourier_ψ::Boson twopi_fourier_χ::Boson

function twopi_fourier_hs_interaction()
    ψc = twopi_fourier_ψ[Classical]
    ψq = twopi_fourier_ψ[Quantum]
    χc = twopi_fourier_χ[Classical]
    χq = twopi_fourier_χ[Quantum]
    forward = ψc^2 * bar(χq) + 2 * ψc * ψq * bar(χc) + (1 // 4) * ψq^2 * bar(χq)
    return ChargedInteractionLagrangian(
        -im * (forward + bar(forward)),
        twopi_fourier_ψ => 1,
        twopi_fourier_χ => 2;
        parameter=:h,
    )
end

function twopi_fourier_derivative_interaction()
    ψc = twopi_fourier_ψ[Classical]
    ψq = twopi_fourier_ψ[Quantum]
    χc = twopi_fourier_χ[Classical]
    χq = twopi_fourier_χ[Quantum]
    forward =
        ψc^2 * bar(χq) + 2 * partial(ψc, :x) * ψq * bar(χc) + (1 // 4) * ψq^2 * bar(χq)
    return ChargedInteractionLagrangian(
        -im * (forward + bar(forward)),
        twopi_fourier_ψ => 1,
        twopi_fourier_χ => 2;
        parameter=:d,
    )
end

function sorted_coefficients(collection)
    values = [coefficient for (_, coefficient) in collection]
    return sort!(values; by=value -> (real(value), imag(value)))
end

function sorted_fourier_coefficients(collection)
    values = [
        contribution.coefficient for (_, contributions) in collection for
        contribution in contributions
    ]
    return sort!(values; by=value -> (real(value), imag(value)))
end

function certify_twopi_fourier_component(coordinate, fourier)
    @test length(coordinate) == length(fourier)
    @test sorted_coefficients(coordinate) == sorted_fourier_coefficients(fourier)
    for (graph, _) in fourier
        @test graph.external_count == 1
        @test graph.loop_count == 1
    end
end

function has_external_momentum_factor(collection)
    for (graph, contributions) in collection
        graph.external_count == 1 || continue
        for contribution in contributions
            for (monomial, _) in contribution.kinematic.terms
                for factor in monomial
                    iszero(factor.momentum[1]) || return true
                end
            end
        end
    end
    return false
end

@testset "pre-amputation 2PI Fourier lowering" begin
    L = twopi_fourier_hs_interaction()
    Γ2 = @inferred TwoPIEffectiveAction(L, Val(2), Val(3))

    for target in (twopi_fourier_ψ, twopi_fourier_χ)
        coordinate = @inferred SelfEnergy(Γ2, target)
        fourier = @inferred KC.twopi_fourier_self_energy(Γ2, target)

        @test KC.parameters(fourier) == KC.parameters(Γ2)
        @test KC.target_family(fourier) == target
        @test KC.order(fourier) == KC.order(Γ2)
        @test KC.statistics(fourier) === Boson

        certify_twopi_fourier_component(
            KC.keldysh_component(coordinate), KC.keldysh_component(fourier)
        )
        certify_twopi_fourier_component(
            KC.retarded_component(coordinate), KC.retarded_component(fourier)
        )
        certify_twopi_fourier_component(
            KC.advanced_component(coordinate), KC.advanced_component(fourier)
        )

        @test_throws ArgumentError fourier_transform(coordinate)
        @test wigner_transform(fourier; gradient_order=Val(0)) !== nothing
    end
end

@testset "2PI cut retains derivative endpoint momentum" begin
    L = twopi_fourier_derivative_interaction()
    Γ2 = @inferred TwoPIEffectiveAction(L, Val(2), Val(3))
    fourier = @inferred KC.twopi_fourier_self_energy(Γ2, twopi_fourier_ψ)

    @test any(
        has_external_momentum_factor,
        (
            KC.keldysh_component(fourier),
            KC.retarded_component(fourier),
            KC.advanced_component(fourier),
        ),
    )
end
