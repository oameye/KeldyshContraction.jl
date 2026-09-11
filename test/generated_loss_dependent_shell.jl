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

function generated_loss_causal_exception_signature(exception)
    term = exception.causal_term
    coincident = KC.coincident_causal_pole_witness(exception)
    zero_energy = KC.zero_energy_causal_denominator_witness(exception)
    return (;
        kind=exception.kind,
        topology=Tuple(KC.topology(exception.state.term)),
        dependent=KC.has_dependent_shell_support(exception),
        active_lines=Tuple(exception.state.active_line_indices),
        active_frequencies=Tuple(exception.state.active_loop_basis_indices),
        state_factor=exception.state.factor,
        support=exception.support,
        causal_coefficient=term.coefficient,
        coincident_frequency=coincident.frequency_index,
        pivot=coincident.pivot_denominator,
        coincident=coincident.coincident_denominator,
        zero_energy_denominator=zero_energy.denominator_index,
        denominators=Tuple(
            (
                Tuple(denominator.loop_coefficients),
                denominator.energy,
                denominator.infinitesimal,
            ) for denominator in term.denominators
        ),
    )
end

@testset "generated γ² dependent shell and Trotter structure" begin
    collision = generated_loss_spectral_collision()
    @test parameters(collision) == KC.ParameterMonomial(:γ)^2
    @test target_family(collision) === dependent_loss_ϕ

    records = NamedTuple[]
    unshifted_dependent_reduced = 0
    regular_sunset_results = 0
    shifted_isolated_trotter = 0
    shifted_nonisolated_trotter = 0
    dependent_shifted_isolated = 0
    nondependent_shifted_isolated = 0

    shifted_regular_results = 0
    shifted_dependent_results = 0
    shifted_zero_results = 0
    shifted_trotter_results = 0
    shifted_causal_source_results = 0
    shifted_causal_branches = 0
    coincident_causal_branches = 0
    zero_energy_causal_branches = 0

    for (part, expression) in (
        (:offset, collision_offset(collision)),
        (:distribution, collision_distribution_coefficient(collision)),
    )
        for (term, coefficient) in expression
            record = generated_loss_record(part, term, coefficient)
            push!(records, record)
            shifted = any(shift -> !iszero(shift), record.all_shifts)

            if shifted
                witness = @inferred KC.trotter_frequency_witness(term)
                if KC.has_isolated_trotter_frequency(witness)
                    shifted_isolated_trotter += 1
                    if record.dependent
                        dependent_shifted_isolated += 1
                    else
                        nondependent_shifted_isolated += 1
                    end
                    line = KC.trotter_frequency_line(term, witness)
                    @test statistical_weight(line) === NoStatisticalWeight
                    @test @inferred(KC.trotter_equal_time_factor(term, witness)) isa
                        KC.ComplexRationals
                else
                    shifted_nonisolated_trotter += 1
                end
            end

            result = @inferred KC.reduce_frequency_term(term, dependent_loss_ϕ)

            if shifted
                !isempty(result.regular) && (shifted_regular_results += 1)
                !isempty(result.dependent) && (shifted_dependent_results += 1)
                !isempty(result.trotter) && (shifted_trotter_results += 1)
                isempty(result) && (shifted_zero_results += 1)
                if !isempty(result.causal)
                    shifted_causal_source_results += 1
                    shifted_causal_branches += length(result.causal)
                    for exception in result.causal
                        if exception.kind === KC.CausalCoincidentPole
                            coincident_causal_branches += 1
                            @test KC.has_coincident_causal_pole(
                                KC.coincident_causal_pole_witness(exception)
                            )
                            @test !KC.has_zero_energy_causal_denominator(
                                KC.zero_energy_causal_denominator_witness(exception)
                            )
                        elseif exception.kind === KC.CausalZeroEnergyDenominator
                            zero_energy_causal_branches += 1
                            @test KC.has_zero_energy_causal_denominator(
                                KC.zero_energy_causal_denominator_witness(exception)
                            )
                            @test !KC.has_coincident_causal_pole(
                                KC.coincident_causal_pole_witness(exception)
                            )
                        else
                            error("unknown generated causal exceptional kind")
                        end
                        @info "post-Trotter causal exception" part coefficient signature=generated_loss_causal_exception_signature(
                            exception
                        )
                    end
                end

                for dependent_term in result.dependent
                    @test !KC.requires_trotter_regularisation(dependent_term)
                    @test !isempty(KC.dependent_residual_support(dependent_term))
                end
                @test isempty(result.trotter)
            elseif record.dependent
                @test isempty(result.regular)
                @test length(result.dependent) == 1
                @test isempty(result.causal)
                @test isempty(result.trotter)
                dependent_term = only(result.dependent)
                @test !KC.requires_trotter_regularisation(dependent_term)
                @test !isempty(KC.dependent_residual_support(dependent_term))
                @test KC.dependent_shell_support(dependent_term) ==
                    KC.dependent_shell_support(
                    KC.analyze_spectral_dependencies(term, dependent_loss_ϕ)
                )
                unshifted_dependent_reduced += 1
            elseif record.topology == (3,) && all(==(CollisionSpectral), record.kinds)
                @test !isempty(result.regular)
                @test isempty(result.dependent)
                @test isempty(result.causal)
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

    @info "second-order pre-Trotter census" shifted_isolated_trotter shifted_nonisolated_trotter dependent_shifted_isolated nondependent_shifted_isolated
    @info "second-order post-Trotter outcomes" shifted_regular_results shifted_dependent_results shifted_zero_results shifted_trotter_results shifted_causal_source_results shifted_causal_branches coincident_causal_branches zero_energy_causal_branches

    @test length(records) == 62
    @test length(dependent_records) == 18
    @test length(shifted_records) == 32
    @test length(unshifted_dependent_records) == 6
    @test unshifted_dependent_reduced == 6
    @test regular_sunset_results == 8

    @test shifted_isolated_trotter == 32
    @test shifted_nonisolated_trotter == 0
    @test dependent_shifted_isolated == 12
    @test nondependent_shifted_isolated == 20
    @test shifted_trotter_results == 0
    @test shifted_causal_branches ==
        coincident_causal_branches + zero_energy_causal_branches
    @test shifted_causal_branches > 0

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
