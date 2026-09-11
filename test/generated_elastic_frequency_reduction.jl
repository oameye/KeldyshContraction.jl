using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields generated_elastic_ϕ::Boson

function generated_elastic_spectral_collision()
    c = generated_elastic_ϕ[Classical]
    q = generated_elastic_ϕ[Quantum]
    interaction =
        -(1 // 2) * ((c^2 + q^2) * bar(c) * bar(q) + c * q * (bar(c)^2 + bar(q)^2))
    L = InteractionLagrangian(interaction, :g)
    G = DressedPropagator(L, Val(2), Val(5); simplify=true)
    GF = fourier_transform(G)
    ΣW = wigner_transform(SelfEnergy(GF); gradient_order=Val(0))
    return spectral_dispersive_collision(
        off_shell_collision_expression(kinetic_expression(ΣW))
    )
end

elastic_momentum_signature(momentum) = Tuple(momentum.coefficients)

function elastic_statistical_signature(term)
    signature = [
        elastic_momentum_signature(KC.momentum(line)) for
        line in kinetic_lines(term.carrier) if
        statistical_weight(line) === DistributionWeight
    ]
    sort!(signature)
    return Tuple(signature)
end

function elastic_reduced_coefficient(records, part, signature)
    coefficient = zero(KC.ComplexRationals)
    for record in records
        record.part === part || continue
        record.signature == signature || continue
        coefficient += record.coefficient
    end
    return coefficient
end

function swap_elastic_qr_momentum(signature)
    length(signature) == 3 ||
        throw(DimensionMismatch("elastic momentum signature must use k,q,r"))
    return (signature[1], signature[3], signature[2])
end

function swap_elastic_qr_signature(signature)
    swapped = [swap_elastic_qr_momentum(momentum) for momentum in signature]
    sort!(swapped)
    return Tuple(swapped)
end

function elastic_qr_symmetrized_coefficient(records, part, signature)
    swapped = swap_elastic_qr_signature(signature)
    return (
        elastic_reduced_coefficient(records, part, signature) +
        elastic_reduced_coefficient(records, part, swapped)
    ) / 2
end

@testset "generated elastic g² collision reduction" begin
    collision = generated_elastic_spectral_collision()
    elastic_parameter = KC.ParameterMonomial(:g)^2
    @test parameters(collision) == elastic_parameter
    @test target_family(collision) === generated_elastic_ϕ

    records = NamedTuple[]
    nterms = 0

    for (part, expression) in (
        (:offset, collision_offset(collision)),
        (:distribution, collision_distribution_coefficient(collision)),
    )
        for (term, coefficient) in expression
            nterms += 1
            @test Tuple(KC.topology(term)) == (3,)
            @test loop_frequency_count(term) == 2
            @test all(iszero ∘ regularisation_shift, kinetic_lines(term.carrier))

            reduction = @inferred KC.general_frequency_reduction(term, generated_elastic_ϕ)
            for (support, frequency_coefficient) in reduction
                push!(
                    records,
                    (;
                        part,
                        term,
                        kinds=Tuple(spectral_dispersive_kinds(term)),
                        rank=spectral_frequency_rank(term),
                        signature=elastic_statistical_signature(term),
                        coefficient=coefficient * frequency_coefficient,
                        support,
                    ),
                )
            end
        end
    end

    @test nterms == 11
    @test length(records) == 11

    # The generated elastic sector must actually exercise the causal fallback. Explicit
    # spectral rank alone is insufficient: rank-one terms with two dispersive lines acquire
    # their second shell constraint through the exact D*D causal convolution.
    @test any(
        record -> record.rank == 1 && count(==(CollisionDispersive), record.kinds) == 2,
        records,
    )

    representative = first(records).term
    basis = KC.momentum_basis(representative)
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    p = -k + q + r
    energy(momentum) = KC.EnergyForm(KC.DispersionAtom(generated_elastic_ϕ, momentum))

    expected_shell, shell_factor = KC.energy_shell(
        energy(k) + energy(p) - energy(q) - energy(r)
    )
    @test shell_factor == 1
    expected_support = KC.FrequencySupport(
        [expected_shell], KC.PrincipalValueSupport{Boson}[]
    )
    @test all(record -> record.support == expected_support, records)

    signature(momentums...) =
        Tuple(sort([elastic_momentum_signature(momentum) for momentum in momentums]))
    pqr = signature(p, q, r)
    p_only = signature(p)
    q_only = signature(q)
    r_only = signature(r)
    pq = signature(p, q)
    pr = signature(p, r)
    qr = signature(q, r)
    none = ()

    expected_signatures = Set((pqr, p_only, q_only, r_only, pq, pr, qr, none))
    @test all(record -> record.signature in expected_signatures, records)

    # The generated routing keeps the two dummy outgoing loop variables q and r distinct. The
    # exact raw affine polynomial is therefore asymmetric before the later loop-momentum
    # quotient:
    #
    #   offset       =  1/2 F_p F_q F_r + 1/2 F_p - F_q
    #   distribution = -F_p F_q + 1/2 F_q F_r + 1/2.
    #
    # These assertions certify what the generated collision actually contains, without
    # silently identifying dummy integration variables.
    @test elastic_reduced_coefficient(records, :offset, pqr) == 1 // 2
    @test elastic_reduced_coefficient(records, :offset, p_only) == 1 // 2
    @test elastic_reduced_coefficient(records, :offset, q_only) == -1
    @test elastic_reduced_coefficient(records, :offset, r_only) == 0
    @test elastic_reduced_coefficient(records, :distribution, pq) == -1
    @test elastic_reduced_coefficient(records, :distribution, pr) == 0
    @test elastic_reduced_coefficient(records, :distribution, qr) == 1 // 2
    @test elastic_reduced_coefficient(records, :distribution, none) == 1 // 2

    # Quotienting only by the legitimate dummy-variable exchange q↔r gives
    #
    # P(F_k,F_p,F_q,F_r)
    #   = 1/2 (F_p F_q F_r + F_p - F_q - F_r)
    #     + F_k/2 (-F_p F_q - F_p F_r + F_q F_r + 1).
    #
    # With F=1+2n this is +4 times the standard Bose gain-loss bracket
    #
    #   (1+n_k)(1+n_p)n_q n_r - n_k n_p(1+n_q)(1+n_r).
    #
    # In the package convention I_pkg = iΣᴷ - iF_k(Σᴿ-Σᴬ), the physical occupation
    # collision is C_n = I_pkg/2. Hence this maps to the standard +2g² Boltzmann kernel.
    @test elastic_qr_symmetrized_coefficient(records, :offset, pqr) == 1 // 2
    @test elastic_qr_symmetrized_coefficient(records, :offset, p_only) == 1 // 2
    @test elastic_qr_symmetrized_coefficient(records, :offset, q_only) == -1 // 2
    @test elastic_qr_symmetrized_coefficient(records, :offset, r_only) == -1 // 2
    @test elastic_qr_symmetrized_coefficient(records, :distribution, pq) == -1 // 2
    @test elastic_qr_symmetrized_coefficient(records, :distribution, pr) == -1 // 2
    @test elastic_qr_symmetrized_coefficient(records, :distribution, qr) == 1 // 2
    @test elastic_qr_symmetrized_coefficient(records, :distribution, none) == 1 // 2
end
