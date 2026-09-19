using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields composite_parity_ψ::Boson composite_parity_χ::Boson

const CompositeParityCoeff = Complex{Rational{Int}}

function composite_parity_hs_interaction()
    ψc = composite_parity_ψ[Classical]
    ψq = composite_parity_ψ[Quantum]
    χc = composite_parity_χ[Classical]
    χq = composite_parity_χ[Quantum]

    forward = (1 // 2) * ψc^2 * bar(χq) + ψc * ψq * bar(χc) + (1 // 2) * ψq^2 * bar(χq)
    return ChargedInteractionLagrangian(
        -im * (forward + bar(forward)),
        composite_parity_ψ => 1,
        composite_parity_χ => 2;
        parameter=:h,
    )
end

function composite_parity_microscopic_interaction()
    c = composite_parity_ψ[Classical]
    q = composite_parity_ψ[Quantum]
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

    return InteractionLagrangian(elastic, :g) + InteractionLagrangian(loss, :γ)
end

function composite_parity_microscopic_skeleton(parameter)
    L = composite_parity_microscopic_interaction()
    G = DressedPropagator(L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false)
    return KC.skeleton_self_energy(fourier_transform(G[parameter]))
end

function composite_parity_component(object, kind::Symbol)
    kind === :K && return KC.keldysh_component(object)
    kind === :R && return KC.retarded_component(object)
    kind === :A && return KC.advanced_component(object)
    return error("unexpected Keldysh component")
end

function composite_parity_component_kind(component::UInt8)
    component == 0x01 && return :R
    component == 0x02 && return :A
    component == 0x03 && return :K
    return error("unexpected physical 2PI component")
end

function composite_parity_line_kind(edge)
    KC.is_keldysh(edge) && return :K
    KC.is_retarded(edge) && return :R
    KC.is_advanced(edge) && return :A
    return error("unexpected bosonic propagator type")
end

function composite_parity_response_kind(kind::KC.PropagatorType.T)
    kind === KC.PropagatorType.Keldysh && return :K
    kind === KC.PropagatorType.Retarded && return :R
    kind === KC.PropagatorType.Advanced && return :A
    return error("unexpected response component")
end

function composite_parity_twopi_cut_terms(Γ2, target)
    records = NamedTuple[]
    derivative_factor = convert(CompositeParityCoeff, im)
    for (vacuum, vacuum_coefficient) in Γ2
        contractions = KC.twopi_contractions(vacuum)
        for cut_index in eachindex(contractions)
            cut = contractions[cut_index]
            KC.field_family(cut.out) == target || continue

            component = KC._twopi_physical_component(
                KC.keldysh_index(cut.in), KC.keldysh_index(cut.out)
            )
            iszero(component) && continue
            internal = KC.Contraction{Boson}[
                contractions[index] for
                index in eachindex(contractions) if index != cut_index
            ]
            KC._twopi_internal_is_physical(internal) || continue

            push!(
                records,
                (;
                    cut,
                    internal,
                    kind=composite_parity_component_kind(component),
                    coefficient=derivative_factor *
                                convert(CompositeParityCoeff, vacuum_coefficient),
                ),
            )
        end
    end
    return records
end

function composite_parity_outer_carrier_graph(outer)
    component = KC._twopi_physical_component(
        KC.keldysh_index(outer.cut.in), KC.keldysh_index(outer.cut.out)
    )
    out_type, in_type = KC._twopi_fourier_external_types(component)
    out_anchor = KC._twopi_fourier_out_anchor(outer.cut, composite_parity_ψ, out_type)
    in_anchor = KC._twopi_fourier_in_anchor(outer.cut, composite_parity_ψ, in_type)
    source = KC.Diagram(
        KC.Contraction{Boson}[outer.internal..., out_anchor, in_anchor], Val(4), Val(1)
    )
    routed, kinematic, _ = KC._canonical_fourier_source(source)
    amputated = KC._amputate_fourier_graph(routed, Val(2), Val(1))
    response_index = KC._composite_response_edge_index(amputated, composite_parity_χ)
    return KC.CompositeFourierDiagram(amputated, response_index), kinematic
end

function composite_parity_expected_carrier_component(Γ2, kind::Symbol)
    result = KC.CompositeFourierDiagrams{CompositeParityCoeff,Boson,2,1}()
    for outer in composite_parity_twopi_cut_terms(Γ2, composite_parity_ψ)
        outer.kind === kind || continue
        graph, kinematic = composite_parity_outer_carrier_graph(outer)
        KC._push_composite_fourier!(result, graph, outer.coefficient, kinematic)
    end
    return result
end

function composite_parity_response_terms(kind::Symbol, sector::Symbol)
    one_c = one(CompositeParityCoeff)
    i_c = complex(0 // 1, 1 // 1)

    if kind === :R
        coefficient = if sector === :g2
            -one_c
        elseif sector === :gγ
            2 * i_c
        elseif sector === :γ2
            one_c
        else
            return error("unexpected weak response sector")
        end
        return ((:R, coefficient),)
    elseif kind === :A
        coefficient = if sector === :g2
            -one_c
        elseif sector === :gγ
            -2 * i_c
        elseif sector === :γ2
            one_c
        else
            return error("unexpected weak response sector")
        end
        return ((:A, coefficient),)
    elseif kind === :K
        sector === :g2 && return ((:K, -one_c),)
        sector === :gγ && return ((:R, 2 * i_c), (:A, 2 * i_c))
        sector === :γ2 && return ((:K, -one_c), (:R, 2 * one_c), (:A, -2 * one_c))
    end
    return error("unexpected weak response component or sector")
end

function composite_parity_reposition(field, row_from, row_to, column_from, column_to)
    p = KC.position(field)
    p == row_from && return field(row_to)
    p == column_from && return field(column_to)
    return error("pair bubble contains a vertex unrelated to its χ cut")
end

function composite_parity_attach_bubble(outer_hs, bubble)
    row_from = KC.position(bubble.cut.in)
    column_from = KC.position(bubble.cut.out)
    row_to = KC.position(outer_hs.out)
    column_to = KC.position(outer_hs.in)

    attached = KC.Contraction{Boson}[]
    sizehint!(attached, length(bubble.internal))
    for contraction in bubble.internal
        push!(
            attached,
            KC.Contraction(
                composite_parity_reposition(
                    contraction.out, row_from, row_to, column_from, column_to
                ),
                composite_parity_reposition(
                    contraction.in, row_from, row_to, column_from, column_to
                ),
            ),
        )
    end
    return attached
end

function composite_parity_anchored_fourier_sunset(outer, contractions)
    length(contractions) == 3 || error("invalid composite weak-expansion sunset")
    component = KC._twopi_physical_component(
        KC.keldysh_index(outer.cut.in), KC.keldysh_index(outer.cut.out)
    )
    out_type, in_type = KC._twopi_fourier_external_types(component)
    out_anchor = KC._twopi_fourier_out_anchor(outer.cut, composite_parity_ψ, out_type)
    in_anchor = KC._twopi_fourier_in_anchor(outer.cut, composite_parity_ψ, in_type)
    source = KC.Diagram(
        KC.Contraction{Boson}[contractions..., out_anchor, in_anchor], Val(5), Val(1)
    )
    routed, kinematic, _ = KC._canonical_fourier_source(source)
    amputated = KC._amputate_fourier_graph(routed, Val(3), Val(1))
    return amputated, kinematic
end

function composite_parity_expanded_component(
    Γ2, carrier, component_kind::Symbol, sector::Symbol
)
    result = KC.FourierDiagrams{CompositeParityCoeff,Boson,3,1}()
    outer_terms = composite_parity_twopi_cut_terms(Γ2, composite_parity_ψ)
    bubble_terms = composite_parity_twopi_cut_terms(Γ2, composite_parity_χ)
    carrier_component = composite_parity_component(carrier, component_kind)

    for outer in outer_terms
        outer.kind === component_kind || continue
        carrier_graph, carrier_kinematic = composite_parity_outer_carrier_graph(outer)
        contributions = get(carrier_component.diagrams, carrier_graph, nothing)
        contributions === nothing &&
            error("production carrier is missing an anchored outer cut")
        any(
            contribution ->
                isequal(contribution.kinematic, carrier_kinematic) &&
                contribution.coefficient == outer.coefficient,
            contributions,
        ) || error("production carrier changed an anchored outer contribution")

        atomic = only(
            contraction for contraction in outer.internal if
            KC.field_family(contraction.out) == composite_parity_ψ
        )
        hs = only(
            contraction for contraction in outer.internal if
            KC.field_family(contraction.out) == composite_parity_χ
        )

        response_kind = composite_parity_response_kind(KC.response_component(carrier_graph))
        response_kind === composite_parity_line_kind(hs) ||
            error("production carrier changed the response Keldysh component")

        for (Ωkind, response_coefficient) in
            composite_parity_response_terms(response_kind, sector)
            for bubble in bubble_terms
                bubble.kind === Ωkind || continue
                contractions = KC.Contraction{Boson}[atomic]
                append!(contractions, composite_parity_attach_bubble(hs, bubble))
                KC._twopi_internal_is_physical(contractions) || continue
                sunset, kinematic = composite_parity_anchored_fourier_sunset(
                    outer, contractions
                )

                coefficient = -outer.coefficient * response_coefficient * bubble.coefficient
                KC._push_fourier!(result, sunset, coefficient, kinematic)
            end
        end
    end
    return result
end

function composite_parity_expanded_self_energy(Γ2, carrier, sector::Symbol, parameter)
    return KC.FourierSelfEnergy{CompositeParityCoeff,Boson,2,3,1}(
        composite_parity_expanded_component(Γ2, carrier, :K, sector),
        composite_parity_expanded_component(Γ2, carrier, :R, sector),
        composite_parity_expanded_component(Γ2, carrier, :A, sector),
        parameter,
        composite_parity_ψ,
    )
end

function composite_parity_canonical_collision(ΣF)
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    KΣ = kinetic_expression(ΣW)
    off_shell = off_shell_collision_expression(KΣ)
    spectral = spectral_dispersive_collision(off_shell)
    return KC.canonical_frequency_collision(spectral)
end

function composite_parity_statistical_families(collision)
    families = Any[]
    for sector in keys(KC.canonical_frequency_expressions(collision))
        for atom in KC.statistical_monomial(sector)
            push!(families, KC.statistical_family(atom))
        end
    end
    for term in KC.shifted_frequency_terms(collision)
        for atom in KC.statistical_monomial(term)
            push!(families, KC.statistical_family(atom))
        end
    end
    return families
end

function composite_parity_support_signature(supports)
    return Tuple(
        (
            KC.affine_support_rank(support),
            KC.affine_support_count(support),
            Tuple(Tuple(row) for row in support.dependency_rows),
        ) for support in supports
    )
end

@testset "production composite carrier equals anchored raw ψ cuts" begin
    Γ2 = TwoPIEffectiveAction(composite_parity_hs_interaction(), Val(2), Val(3))
    carrier = KC.twopi_composite_fourier_self_energy(
        Γ2, composite_parity_ψ, composite_parity_χ
    )

    for kind in (:K, :R, :A)
        @test isequal(
            composite_parity_component(carrier, kind),
            composite_parity_expected_carrier_component(Γ2, kind),
        )
    end
end

@testset "production composite carrier reproduces frozen C3 weak parity" begin
    Γ2 = TwoPIEffectiveAction(composite_parity_hs_interaction(), Val(2), Val(3))
    carrier = KC.twopi_composite_fourier_self_energy(
        Γ2, composite_parity_ψ, composite_parity_χ
    )
    sectors = (
        (:g2, KC.ParameterMonomial(:g)^2),
        (:gγ, KC.ParameterMonomial(:g) * KC.ParameterMonomial(:γ)),
        (:γ2, KC.ParameterMonomial(:γ)^2),
    )

    expanded = Dict{Symbol,Any}()
    direct = Dict{Symbol,Any}()
    for (sector, parameter) in sectors
        expanded[sector] = composite_parity_expanded_self_energy(
            Γ2, carrier, sector, parameter
        )
        direct[sector] = composite_parity_microscopic_skeleton(parameter)

        for component_kind in (:K, :R, :A)
            component = composite_parity_component(expanded[sector], component_kind)
            @test !isempty(component)
            for (graph, contributions) in component
                @test graph.external_count == 1
                @test graph.loop_count == 2
                @test collect(KC.topology(graph.coordinate)) == [3]
                @test length(KC.contractions(graph.coordinate)) == 3
                @test all(
                    edge -> KC.field_family(edge.out) == composite_parity_ψ,
                    KC.contractions(graph.coordinate),
                )
                @test !isempty(contributions)
            end
        end
    end

    @test isequal(expanded[:g2].keldysh, direct[:g2].keldysh)
    @test isequal(expanded[:g2].retarded, direct[:g2].retarded)
    @test isequal(expanded[:g2].advanced, direct[:g2].advanced)

    mixed_expanded = composite_parity_canonical_collision(expanded[:gγ])
    mixed_direct = composite_parity_canonical_collision(direct[:gγ])
    @test isempty(KC.shifted_frequency_terms(mixed_expanded))
    @test isempty(KC.shifted_frequency_terms(mixed_direct))
    @test KC.canonical_frequency_expressions(mixed_expanded) ==
        KC.canonical_frequency_expressions(mixed_direct)
    @test all(==(composite_parity_ψ), composite_parity_statistical_families(mixed_expanded))

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

    expanded_blockers = KC.blocked_trotter_contributions(expanded_reduction)
    direct_blockers = KC.blocked_trotter_contributions(direct_reduction)
    @test Set(keys(expanded_blockers)) == Set(keys(direct_blockers))
    for (key, expanded_blocked) in expanded_blockers
        direct_blocked = direct_blockers[key]
        @test KC.trotter_frequency_blocker_kind(expanded_blocked) ===
            KC.trotter_frequency_blocker_kind(direct_blocked)
        @test KC.trotter_frequency_blocked_index(expanded_blocked) ==
            KC.trotter_frequency_blocked_index(direct_blocked)
        @test KC.trotter_frequency_blocked_expression(expanded_blocked) ==
            KC.trotter_frequency_blocked_expression(direct_blocked)
        @test composite_parity_support_signature(
            KC.trotter_frequency_blocked_supports(expanded_blocked)
        ) == composite_parity_support_signature(
            KC.trotter_frequency_blocked_supports(direct_blocked)
        )
    end
end
