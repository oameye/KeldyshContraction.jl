using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields dependent_loss_ϕ::Boson

function generated_loss_spectral_collision()
    c = dependent_loss_ϕ[Classical]
    q = dependent_loss_ϕ[Quantum]
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
    G = DressedPropagator(L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false)
    GF = fourier_transform(G)
    ΣW = wigner_transform(SelfEnergy(GF); gradient_order=Val(0))
    collision = off_shell_collision_expression(kinetic_expression(ΣW))
    return spectral_dispersive_collision(collision)
end

function generated_loss_record(part, term, coefficient)
    analysis = @inferred KC.analyze_spectral_dependencies(term, dependent_loss_ϕ)
    support = KC.dependent_shell_support(analysis)
    all_shifts = Tuple(regularisation_shift(line) for line in kinetic_lines(term.carrier))
    return (;
        part,
        topology=Tuple(KC.topology(term)),
        kinds=Tuple(spectral_dispersive_kinds(term)),
        coefficient,
        rank=KC.constraint_rank(analysis),
        count=KC.constraint_count(analysis),
        dependent=KC.has_dependent_shell_support(analysis),
        multiplicity=KC.repeated_shell_multiplicity(analysis),
        dependencies=Tuple(Tuple(row) for row in support.dependency_rows),
        spectral_shifts=Tuple(analysis.regularisation_shifts),
        all_shifts,
    )
end

@testset "generated γ² dependent shell structure" begin
    collision = generated_loss_spectral_collision()
    @test parameters(collision) == KC.ParameterMonomial(:γ)^2
    @test target_family(collision) === dependent_loss_ϕ

    records = NamedTuple[]
    dependent_results = 0
    dependent_reduced_results = 0
    trotter_results = 0
    regular_sunset_results = 0

    for (part, expression) in (
        (:offset, collision_offset(collision)),
        (:distribution, collision_distribution_coefficient(collision)),
    )
        for (term, coefficient) in expression
            record = generated_loss_record(part, term, coefficient)
            push!(records, record)

            if record.dependent
                result = @inferred KC.reduce_frequency_term(term, dependent_loss_ϕ)
                @test isempty(result.regular)
                @test length(result.dependent) == 1
                @test isempty(result.trotter)
                dependent_term = first(result.dependent)
                @test KC.dependent_shell_support(dependent_term) ==
                    KC.dependent_shell_support(
                    KC.analyze_spectral_dependencies(term, dependent_loss_ϕ)
                )
                if any(shift -> !iszero(shift), record.all_shifts)
                    @test KC.requires_trotter_regularisation(dependent_term)
                    @test isempty(KC.dependent_residual_support(dependent_term))
                else
                    @test !KC.requires_trotter_regularisation(dependent_term)
                    @test !isempty(KC.dependent_residual_support(dependent_term))
                    dependent_reduced_results += 1
                end
                dependent_results += 1
            elseif any(shift -> !iszero(shift), record.all_shifts)
                result = @inferred KC.reduce_frequency_term(term, dependent_loss_ϕ)
                @test isempty(result.regular)
                @test isempty(result.dependent)
                @test length(result.trotter) == 1
                trotter_results += 1
            elseif record.topology == (3,) && all(==(CollisionSpectral), record.kinds)
                result = @inferred KC.reduce_frequency_term(term, dependent_loss_ϕ)
                @test !isempty(result.regular)
                @test isempty(result.dependent)
                @test isempty(result.trotter)
                regular_sunset_results += 1
            end
        end
    end

    dependent_records = filter(record -> record.dependent, records)
    shifted_records = filter(
        record -> any(shift -> !iszero(shift), record.all_shifts), records
    )
    unshifted_dependent_records = filter(
        record -> record.dependent && all(iszero, record.all_shifts), records
    )

    @test length(records) == 62
    @test length(dependent_records) == 18
    @test length(shifted_records) == 32
    @test length(unshifted_dependent_records) == 6
    @test dependent_results == 18
    @test dependent_reduced_results == 6
    @test trotter_results == 20
    @test regular_sunset_results == 8

    @test all(record -> record.topology == (2,), dependent_records)
    @test all(record -> record.topology == (2,), unshifted_dependent_records)
    @test all(record -> record.rank == 2, unshifted_dependent_records)
    @test all(record -> record.count == 3, unshifted_dependent_records)
    @test all(record -> record.multiplicity == 0, unshifted_dependent_records)
    @test all(
        record -> record.dependencies == ((0 // 1, 1 // 1),), unshifted_dependent_records
    )
    @test all(
        record -> all(==(CollisionSpectral), record.kinds), unshifted_dependent_records
    )
    @test any(record -> record.topology == (2,), shifted_records)
    @test any(record -> record.topology == (3,), records)
end
