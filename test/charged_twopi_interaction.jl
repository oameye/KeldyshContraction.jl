using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields charged_twopi_ψ::Boson charged_twopi_χ::Boson

function charged_twopi_hs_interaction()
    ψc = charged_twopi_ψ[Classical]
    ψq = charged_twopi_ψ[Quantum]
    χq = charged_twopi_χ[Quantum]

    forward = ψc^2 * bar(χq) + 2 * ψc * ψq * bar(χq) + (1 // 4) * ψq^2 * bar(χq)
    interaction = -im * (forward + bar(forward))
    return ChargedInteractionLagrangian(
        interaction, charged_twopi_ψ => 1, charged_twopi_χ => 2; parameter=:h
    )
end

function charged_twopi_is_line(contraction, family, component)
    return KC.field_family(contraction.out) == family &&
           KC.keldysh_index(contraction.out) === first(component) &&
           KC.keldysh_index(contraction.in) === last(component)
end

@testset "charged interaction preserves weighted U(1) neutrality" begin
    L = charged_twopi_hs_interaction()

    @test KC.field_families(L) == sort([charged_twopi_ψ, charged_twopi_χ])
    @test KC.field_charge(L, charged_twopi_ψ) == 1
    @test KC.field_charge(L, charged_twopi_χ) == 2
    @test KC.parameters(L) == KC.parameter_monomial(:h)
    @test_throws ArgumentError KC.target_family(L)
    @test KC.target_family(L, charged_twopi_ψ) == charged_twopi_ψ

    ψc = charged_twopi_ψ[Classical]
    ψq = charged_twopi_ψ[Quantum]
    χq = charged_twopi_χ[Quantum]
    forward = ψc^2 * bar(χq) + 2 * ψc * ψq * bar(χq) + (1 // 4) * ψq^2 * bar(χq)
    interaction = -im * (forward + bar(forward))

    # The existing interaction contract remains deliberately stricter.
    @test_throws ArgumentError InteractionLagrangian(interaction, :h)
    @test_throws ArgumentError ChargedInteractionLagrangian(
        interaction, charged_twopi_ψ => 1, charged_twopi_χ => 1; parameter=:h
    )
end

@testset "charged cubic interaction generates the G²D 2PI skeleton" begin
    L = charged_twopi_hs_interaction()
    Γ2 = @inferred TwoPIEffectiveAction(L, Val(2), Val(3))

    @test !iszero(Γ2)
    @test KC.order(Γ2) == 2
    @test KC.parameters(Γ2) == KC.parameters(L)^2
    @test KC.field_families(Γ2) == KC.field_families(L)

    terms = KC.twopi_terms(Γ2)
    @test !isempty(terms)
    @test length(unique(KC.twopi_topology(diagram) for diagram in keys(terms))) == 1

    for diagram in keys(terms)
        families = [
            KC.field_family(contraction.out) for
            contraction in KC.twopi_contractions(diagram)
        ]
        @test count(==(charged_twopi_ψ), families) == 2
        @test count(==(charged_twopi_χ), families) == 1
        @test KC.is_two_particle_irreducible(collect(KC.twopi_contractions(diagram)))
    end

    # This component comes only from ψ_c² χ̄_q times its conjugate. KC's perturbative
    # and propagator phase conventions give the exact stored coefficient -2 h²: the
    # factor two is the two equivalent ψ Wick pairings. The source's continuum
    # -i h² G²D prefactor is recovered only after translating propagator conventions.
    ccqq_terms = [
        (diagram, coefficient) for (diagram, coefficient) in terms if begin
            contractions = KC.twopi_contractions(diagram)
            count(
                contraction -> charged_twopi_is_line(
                    contraction, charged_twopi_ψ, (Classical, Classical)
                ),
                contractions,
            ) == 2 &&
                count(
                    contraction -> charged_twopi_is_line(
                        contraction, charged_twopi_χ, (Quantum, Quantum)
                    ),
                    contractions,
                ) == 1
        end
    ]
    @test length(ccqq_terms) == 1
    @test last(only(ccqq_terms)) == -2

    Σψ = @inferred KC.twopi_self_energy(Γ2, charged_twopi_ψ)
    Ωχ = @inferred KC.twopi_self_energy(Γ2, charged_twopi_χ)
    @test !iszero(Σψ)
    @test !iszero(Ωχ)

    for diagram in keys(KC.twopi_self_energy_terms(Σψ))
        families = [
            KC.field_family(contraction.out) for
            contraction in KC.twopi_self_energy_contractions(diagram)
        ]
        @test count(==(charged_twopi_ψ), families) == 1
        @test count(==(charged_twopi_χ), families) == 1
    end

    for diagram in keys(KC.twopi_self_energy_terms(Ωχ))
        families = [
            KC.field_family(contraction.out) for
            contraction in KC.twopi_self_energy_contractions(diagram)
        ]
        @test count(==(charged_twopi_ψ), families) == 2
        @test count(==(charged_twopi_χ), families) == 0
    end

    # The same graph-cut engine must reproduce the 2:1 functional multiplicity of
    # G²D: cutting either identical ψ line gives Σ ~ 2GD, while cutting the single
    # χ line gives Ω ~ G². These are formal coefficients before physical projection.
    ψ_cc_terms = [
        coefficient for (diagram, coefficient) in KC.twopi_self_energy_terms(Σψ) if
        KC.twopi_self_energy_component(diagram) == (Classical, Classical) && begin
            contractions = KC.twopi_self_energy_contractions(diagram)
            count(
                contraction -> charged_twopi_is_line(
                    contraction, charged_twopi_ψ, (Classical, Classical)
                ),
                contractions,
            ) == 1 &&
                count(
                    contraction -> charged_twopi_is_line(
                        contraction, charged_twopi_χ, (Quantum, Quantum)
                    ),
                    contractions,
                ) == 1
        end
    ]
    χ_qq_terms = [
        coefficient for (diagram, coefficient) in KC.twopi_self_energy_terms(Ωχ) if
        KC.twopi_self_energy_component(diagram) == (Quantum, Quantum) && all(
            contraction -> charged_twopi_is_line(
                contraction, charged_twopi_ψ, (Classical, Classical)
            ),
            KC.twopi_self_energy_contractions(diagram),
        )
    ]
    @test ψ_cc_terms == [-4im]
    @test χ_qq_terms == [-2im]
    @test only(ψ_cc_terms) == 2 * only(χ_qq_terms)
end
