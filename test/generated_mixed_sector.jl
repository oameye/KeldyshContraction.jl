using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields mixed_ϕ::Boson

function generated_mixed_off_shell_collision()
    c = mixed_ϕ[Classical]
    q = mixed_ϕ[Quantum]
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus

    elastic = -(1 // 2) * ((c^2 + q^2) * bar(c) * bar(q) + c * q * (bar(c)^2 + bar(q)^2))
    loss =
        (1 // 2) *
        im *
        (
            bar(c) * bar(q) * (c(minus)^2 + q(minus)^2) -
            c(plus) * q(plus) * (bar(c)^2 + bar(q)^2) +
            2 * bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
        )

    L = InteractionLagrangian(elastic, :g) + InteractionLagrangian(loss, :γ)
    G = DressedPropagator(L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false)
    mixed_parameter = KC.ParameterMonomial(:g) * KC.ParameterMonomial(:γ)
    GF = fourier_transform(G[mixed_parameter])
    ΣF = SelfEnergy(GF)
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    KΣ = kinetic_expression(ΣW)
    return off_shell_collision_expression(KΣ)
end

momentum_signature(momentum) = Tuple(momentum.coefficients)

function statistical_momentum_signature(lines)
    signature = [
        momentum_signature(KC.momentum(line)) for
        line in lines if statistical_weight(line) === DistributionWeight
    ]
    sort!(signature)
    return Tuple(signature)
end

function reduced_coefficient(records, part, signature)
    coefficient = zero(KC.ComplexRationals)
    for record in records
        record.part === part || continue
        record.signature == signature || continue
        coefficient += record.coefficient
    end
    return coefficient
end

@testset "generated mixed gγ collision reduction" begin
    off_shell = generated_mixed_off_shell_collision()
    mixed_parameter = KC.ParameterMonomial(:g) * KC.ParameterMonomial(:γ)
    @test parameters(off_shell) == mixed_parameter
    @test target_family(off_shell) === mixed_ϕ

    # Before R/A decomposition, topology [2] is already the isolated repeated-frequency
    # causal structure. In each affine collision component there are exactly two terms, one
    # retarded and one advanced, with the same coefficient and the same physical F weights.
    # The causal line shares its routed momentum with one spectral line; the third line carries
    # the independent tadpole momentum. Thus the q⁰ integral factors as
    # A_q (Gᴿ_q + Gᴬ_q), whose causal principal-value part is the Kramers--Kronig zero.
    rak_lollipop = NamedTuple[]
    for (part, expression) in (
        (:offset, collision_offset(off_shell)),
        (:distribution, collision_distribution_coefficient(off_shell)),
    )
        for (term, coefficient) in expression
            collect(term.topology) == [2] || continue
            lines = kinetic_lines(term)
            @test all(iszero ∘ regularisation_shift, lines)

            causal = [
                line for line in lines if
                kinetic_line_kind(line) in (KineticRetarded, KineticAdvanced)
            ]
            spectral = [
                line for line in lines if kinetic_line_kind(line) === KineticSpectral
            ]
            @test length(causal) == 1
            @test length(spectral) == 2

            causal_line = only(causal)
            repeated = KC.momentum(causal_line)
            @test count(line -> KC.momentum(line) == repeated, spectral) == 1
            @test count(line -> KC.momentum(line) != repeated, spectral) == 1

            push!(
                rak_lollipop,
                (;
                    part,
                    causal_kind=kinetic_line_kind(causal_line),
                    repeated=momentum_signature(repeated),
                    signature=statistical_momentum_signature(lines),
                    coefficient,
                ),
            )
        end
    end

    @test length(rak_lollipop) == 4
    for part in (:offset, :distribution)
        records = filter(record -> record.part === part, rak_lollipop)
        expected_coefficient = part === :offset ? 2 : -2
        @test length(records) == 2
        @test Set(record.causal_kind for record in records) ==
            Set((KineticRetarded, KineticAdvanced))
        @test all(record -> record.coefficient == expected_coefficient, records)
        @test records[1].repeated == records[2].repeated
        @test records[1].signature == records[2].signature
    end

    collision = @inferred spectral_dispersive_collision(off_shell)
    @test parameters(collision) == mixed_parameter
    @test target_family(collision) === mixed_ϕ

    lollipop_kk = 0
    lollipop_spectral = NamedTuple[]
    sunset = NamedTuple[]
    nterms = 0

    for (part, expression) in (
        (:offset, collision_offset(collision)),
        (:distribution, collision_distribution_coefficient(collision)),
    )
        for (term, coefficient) in expression
            nterms += 1
            @test all(iszero ∘ regularisation_shift, kinetic_lines(term.carrier))
            topology = collect(KC.topology(term))
            kinds = spectral_dispersive_kinds(term)

            if topology == [2]
                if count(==(CollisionDispersive), kinds) == 1
                    classification = @inferred classify_exceptional_frequency(term)
                    @test exceptional_frequency_kind(classification) ===
                        FrequencyKramersKronigZero
                    @test exceptional_frequency_loop_basis_index(classification) != 0
                    lollipop_kk += 1
                else
                    @test all(==(CollisionSpectral), kinds)
                    classification = @inferred classify_exceptional_frequency(term)
                    @test exceptional_frequency_kind(classification) === FrequencyUnresolved
                    push!(
                        lollipop_spectral,
                        (;
                            part,
                            signature=statistical_momentum_signature(
                                kinetic_lines(term.carrier)
                            ),
                            coefficient,
                        ),
                    )
                end
            elseif topology == [3]
                reduction = @inferred full_rank_frequency_reduction(term, mixed_ϕ)
                push!(
                    sunset,
                    (;
                        part,
                        term,
                        reduction,
                        signature=statistical_momentum_signature(
                            kinetic_lines(term.carrier)
                        ),
                        coefficient=coefficient * frequency_factor(reduction),
                    ),
                )
            else
                @test false
            end
        end
    end

    @test nterms == 14
    @test lollipop_kk == 4
    @test length(lollipop_spectral) == 4
    @test length(sunset) == 6

    # The all-spectral lollipop pieces differ only in which copy of the repeated q line
    # carries F_q. Once statistical identity is keyed by routed momentum, the two pieces in
    # each affine collision component cancel exactly. Together with the four isolated A_q D_q
    # Kramers--Kronig zeros above, the complete generated topology-[2] sector is zero.
    for part in (:offset, :distribution)
        records = filter(record -> record.part === part, lollipop_spectral)
        @test length(records) == 2
        @test records[1].signature == records[2].signature
        @test sum(record.coefficient for record in records) == 0
    end

    # Every topology-[3] term reduces onto the same four-energy principal-value channel.
    support = frequency_support(first(sunset).reduction)
    @test all(record -> frequency_support(record.reduction) == support, sunset)
    @test isempty(support.shells)
    @test length(support.principal_values) == 1

    representative = first(sunset).term
    basis = KC.momentum_basis(representative)
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    p = -k + q + r
    energy(momentum) = KC.EnergyForm(KC.DispersionAtom(mixed_ϕ, momentum))
    expected_pv, _ = KC.principal_value_support(
        energy(k) + energy(p) - energy(q) - energy(r)
    )
    @test only(support.principal_values) == expected_pv

    q_signature = (momentum_signature(q),)
    qp_signature = Tuple(sort([momentum_signature(q), momentum_signature(p)]))

    # Exact affine F-polynomial in the package collision convention:
    #
    #   offset       =  4 F_q F_p - 4 F_q
    #   distribution = -4 F_q F_p + 4 F_q
    #
    # hence I_gγ = 4 F_q (F_p - 1) (1 - F_k) PV(1/ΔE),
    # ΔE = ε_k + ε_p - ε_q - ε_r and p = -k + q + r.
    # With F=1+2n this is -32 n_k n_p (n_q+1/2) PV(1/ΔE).
    # The independent 2PI convention has C_2PI = -I_package/2, fixed by the elastic g²
    # gain--loss kernel (-4 here versus +2 in 2PI), and therefore gives the independently
    # derived +16 gγ anomalous coefficient.
    @test reduced_coefficient(sunset, :offset, qp_signature) == 4
    @test reduced_coefficient(sunset, :offset, q_signature) == -4
    @test reduced_coefficient(sunset, :distribution, qp_signature) == -4
    @test reduced_coefficient(sunset, :distribution, q_signature) == 4
end
