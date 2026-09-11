using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields generated_loss1_ϕ::Boson

function generated_first_order_loss_collision()
    c = generated_loss1_ϕ[Classical]
    q = generated_loss1_ϕ[Quantum]
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
    G = DressedPropagator(L, Val(1), Val(3); simplify=true, _set_reg_to_zero=false)
    GF = fourier_transform(G)
    ΣW = wigner_transform(SelfEnergy(GF); gradient_order=Val(0))
    return spectral_dispersive_collision(
        off_shell_collision_expression(kinetic_expression(ΣW))
    )
end

loss1_momentum_signature(momentum) = Tuple(momentum.coefficients)

function loss1_statistical_signature(term)
    signature = [
        loss1_momentum_signature(KC.momentum(line)) for
        line in kinetic_lines(term.carrier) if
        statistical_weight(line) === DistributionWeight
    ]
    sort!(signature)
    return Tuple(signature)
end

function loss1_reduced_coefficient(records, part, signature)
    coefficient = zero(KC.ComplexRationals)
    for record in records
        record.part === part || continue
        record.signature == signature || continue
        coefficient += record.coefficient
    end
    return coefficient
end

@testset "generated first-order γ collision reduction" begin
    collision = generated_first_order_loss_collision()
    @test parameters(collision) == KC.ParameterMonomial(:γ)
    @test target_family(collision) === generated_loss1_ϕ

    records = NamedTuple[]
    zero_terms = 0
    source_terms = 0
    shifted_source_terms = 0

    for (part, expression) in (
        (:offset, collision_offset(collision)),
        (:distribution, collision_distribution_coefficient(collision)),
    )
        for (term, coefficient) in expression
            source_terms += 1
            shifted = any(
                line -> !iszero(regularisation_shift(line)), kinetic_lines(term.carrier)
            )
            shifted && (shifted_source_terms += 1)

            result = @inferred KC.reduce_frequency_term(term, generated_loss1_ϕ)
            @test isempty(result.dependent)
            @test isempty(result.causal)
            @test isempty(result.trotter)
            if isempty(result.regular)
                zero_terms += 1
                continue
            end

            signature = loss1_statistical_signature(term)
            for (support, frequency_coefficient) in result.regular
                push!(
                    records,
                    (;
                        part,
                        term,
                        shifted,
                        signature,
                        coefficient=coefficient * frequency_coefficient,
                        support,
                    ),
                )
            end
        end
    end

    grouped = Dict{Tuple{Symbol,Tuple},KC.ComplexRationals}()
    for record in records
        key = (record.part, record.signature)
        grouped[key] = get(grouped, key, zero(KC.ComplexRationals)) + record.coefficient
    end

    @test source_terms == 6
    @test shifted_source_terms == 4
    @test zero_terms == 0
    @test length(grouped) == 4
    @test all(
        record ->
            isempty(record.support.shells) && isempty(record.support.principal_values),
        records,
    )

    # Independent first-order pair-loss oracle:
    #
    #   I_package^γ / γ = 2 (F_k - 1) (1 - F_q)
    #                   = (-2 + 2F_q) + F_k (2 - 2F_q).
    #
    # `part == :offset` stores the F_k-independent term and `:distribution` its F_k
    # coefficient. The internal loop momentum is intentionally inferred from the generated
    # statistical signature rather than hard-coding a basis index.
    internal_signatures = Set(key[2] for key in keys(grouped) if !isempty(key[2]))
    @test length(internal_signatures) == 1
    q_signature = only(internal_signatures)
    @test length(q_signature) == 1

    @test get(grouped, (:offset, ()), zero(KC.ComplexRationals)) == -2 // 1
    @test get(grouped, (:offset, q_signature), zero(KC.ComplexRationals)) == 2 // 1
    @test get(grouped, (:distribution, ()), zero(KC.ComplexRationals)) == 2 // 1
    @test get(grouped, (:distribution, q_signature), zero(KC.ComplexRationals)) == -2 // 1

    # Collision-level assembly must now recover the same oracle without presentation-layer
    # grouping and without using topology as part of the physical term identity.
    reduced = @inferred reduce_frequency_collision(collision)
    @test parameters(reduced) == KC.ParameterMonomial(:γ)
    @test target_family(reduced) === generated_loss1_ϕ
    @test isempty(reduced_dependent_terms(reduced))
    @test isempty(reduced_causal_terms(reduced))
    @test isempty(reduced_trotter_terms(reduced))
    @test length(reduced_regular_terms(reduced)) == 1

    sector, statistical = only(reduced_regular_terms(reduced))
    @test parameters(sector) == KC.ParameterMonomial(:γ)
    @test isempty(frequency_support(sector).shells)
    @test isempty(frequency_support(sector).principal_values)

    basis = KC.momentum_basis(sector)
    external_variable = external_wigner_momentum(sector)
    external_index = only(
        i for (i, variable) in enumerate(basis) if variable == external_variable
    )
    k_momentum = KC.basis_momentum(basis, external_index)
    Fk_atom = StatisticalAtom(generated_loss1_ϕ, k_momentum)

    atoms = Set(atom for (monomial, _) in statistical for atom in monomial)
    @test Fk_atom in atoms
    internal_atoms = [atom for atom in atoms if atom != Fk_atom]
    @test length(internal_atoms) == 1
    Fq_atom = only(internal_atoms)

    Fk = StatisticalPolynomial(Fk_atom, one(KC.ComplexRationals))
    Fq = StatisticalPolynomial(Fq_atom, one(KC.ComplexRationals))
    expected_F = 2 * (Fk - one(typeof(Fk))) * (one(typeof(Fq)) - Fq)
    @test statistical == expected_F

    occupation = @inferred occupation_reduced_expression(reduced)
    @test length(occupation_reduced_terms(occupation)) == 1
    occupation_sector, occupation_polynomial = only(occupation_reduced_terms(occupation))
    @test occupation_sector == sector

    nk = OccupationPolynomial(
        OccupationAtom(generated_loss1_ϕ, k_momentum), one(KC.ComplexRationals)
    )
    nq = OccupationPolynomial(
        OccupationAtom(generated_loss1_ϕ, KC.momentum(Fq_atom)), one(KC.ComplexRationals)
    )
    @test occupation_polynomial == -4 * nk * nq
end
