using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields c5_weak_ψ::Boson c5_weak_χ::Boson

const C5WeakCoeff = Complex{Rational{Int}}

function c5_weak_hs_interaction()
    ψc = c5_weak_ψ[Classical]
    ψq = c5_weak_ψ[Quantum]
    χc = c5_weak_χ[Classical]
    χq = c5_weak_χ[Quantum]
    forward = (1 // 2) * ψc^2 * bar(χq) + ψc * ψq * bar(χc) + (1 // 2) * ψq^2 * bar(χq)
    return ChargedInteractionLagrangian(
        -im * (forward + bar(forward)), c5_weak_ψ => 1, c5_weak_χ => 2; parameter=:h
    )
end

function c5_weak_microscopic_interaction()
    c = c5_weak_ψ[Classical]
    q = c5_weak_ψ[Quantum]
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

function c5_weak_direct_skeleton(parameter)
    G = DressedPropagator(
        c5_weak_microscopic_interaction(),
        Val(2),
        Val(5);
        simplify=true,
        _set_reg_to_zero=false,
    )
    return KC.skeleton_self_energy(fourier_transform(G[parameter]))
end

function c5_weak_component(object, kind::Symbol)
    kind === :K && return KC.keldysh_component(object)
    kind === :R && return KC.retarded_component(object)
    kind === :A && return KC.advanced_component(object)
    return error("unexpected Keldysh component")
end

function c5_weak_component_kind(component::UInt8)
    component == 0x01 && return :R
    component == 0x02 && return :A
    component == 0x03 && return :K
    return error("unexpected physical 2PI component")
end

function c5_weak_response_kind(kind::KC.PropagatorType.T)
    kind === KC.PropagatorType.Keldysh && return :K
    kind === KC.PropagatorType.Retarded && return :R
    kind === KC.PropagatorType.Advanced && return :A
    return error("unexpected response component")
end

function c5_weak_line_kind(edge)
    KC.is_keldysh(edge) && return :K
    KC.is_retarded(edge) && return :R
    KC.is_advanced(edge) && return :A
    return error("unexpected bosonic propagator type")
end

function c5_weak_cut_terms(Γ2, target)
    records = NamedTuple[]
    derivative_factor = convert(C5WeakCoeff, im)
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
                    kind=c5_weak_component_kind(component),
                    coefficient=derivative_factor *
                                convert(C5WeakCoeff, vacuum_coefficient),
                ),
            )
        end
    end
    return records
end

function c5_weak_outer_carrier_graph(outer)
    component = KC._twopi_physical_component(
        KC.keldysh_index(outer.cut.in), KC.keldysh_index(outer.cut.out)
    )
    out_type, in_type = KC._twopi_fourier_external_types(component)
    out_anchor = KC._twopi_fourier_out_anchor(outer.cut, c5_weak_ψ, out_type)
    in_anchor = KC._twopi_fourier_in_anchor(outer.cut, c5_weak_ψ, in_type)
    source = KC.Diagram(
        KC.Contraction{Boson}[outer.internal..., out_anchor, in_anchor], Val(4), Val(1)
    )
    routed, kinematic, _ = KC._canonical_fourier_source(source)
    amputated = KC._amputate_fourier_graph(routed, Val(2), Val(1))
    response_index = KC._composite_response_edge_index(amputated, c5_weak_χ)
    return KC.CompositeFourierDiagram(amputated, response_index), kinematic
end

function c5_weak_response_terms(kind::Symbol, sector::Symbol)
    one_c = one(C5WeakCoeff)
    i_c = complex(0 // 1, 1 // 1)
    kind === :R && return ((
        :ΩR,
        if sector === :g2
            -one_c
        elseif sector === :gγ
            2 * i_c
        else
            one_c
        end,
    ),)
    kind === :A && return ((
        :ΩA,
        if sector === :g2
            -one_c
        elseif sector === :gγ
            -2 * i_c
        else
            one_c
        end,
    ),)
    kind === :K && sector === :g2 && return ((:ΩK, -one_c),)
    kind === :K && sector === :gγ && return ((:ΩR, 2 * i_c), (:ΩA, 2 * i_c))
    kind === :K &&
        sector === :γ2 &&
        return ((:ΩK, -one_c), (:ΩR, 2 * one_c), (:ΩA, -2 * one_c))
    return error("unexpected weak response sector")
end

c5_weak_Ωkind(kind::Symbol) =
    if kind === :ΩK
        :K
    elseif kind === :ΩR
        :R
    elseif kind === :ΩA
        :A
    else
        error("unexpected polarization component")
    end

function c5_weak_reposition(field, row_from, row_to, column_from, column_to)
    position = KC.position(field)
    position == row_from && return field(row_to)
    position == column_from && return field(column_to)
    return error("pair bubble contains an unrelated vertex")
end

function c5_weak_attach_bubble(outer_hs, bubble)
    row_from = KC.position(bubble.cut.in)
    column_from = KC.position(bubble.cut.out)
    row_to = KC.position(outer_hs.out)
    column_to = KC.position(outer_hs.in)
    attached = KC.Contraction{Boson}[]
    for contraction in bubble.internal
        push!(
            attached,
            KC.Contraction(
                c5_weak_reposition(
                    contraction.out, row_from, row_to, column_from, column_to
                ),
                c5_weak_reposition(
                    contraction.in, row_from, row_to, column_from, column_to
                ),
            ),
        )
    end
    return attached
end

function c5_weak_anchored_sunset(outer, contractions)
    component = KC._twopi_physical_component(
        KC.keldysh_index(outer.cut.in), KC.keldysh_index(outer.cut.out)
    )
    out_type, in_type = KC._twopi_fourier_external_types(component)
    out_anchor = KC._twopi_fourier_out_anchor(outer.cut, c5_weak_ψ, out_type)
    in_anchor = KC._twopi_fourier_in_anchor(outer.cut, c5_weak_ψ, in_type)
    source = KC.Diagram(
        KC.Contraction{Boson}[contractions..., out_anchor, in_anchor], Val(5), Val(1)
    )
    routed, kinematic, _ = KC._canonical_fourier_source(source)
    return KC._amputate_fourier_graph(routed, Val(3), Val(1)), kinematic
end

function c5_weak_wigner_graph(component, carrier_graph, carrier_kinematic)
    matches = Tuple[]
    for (graph, contributions) in component
        isequal(graph.composite, carrier_graph) || continue
        any(c -> isequal(c.kinematic, carrier_kinematic), contributions) || continue
        push!(matches, (graph, contributions))
    end
    length(matches) == 1 ||
        error("response-aware Wigner carrier changed anchored outer provenance")
    return only(matches)[1]
end

function c5_weak_expand_component(
    Γ2, response_wigner, component_kind::Symbol, sector::Symbol
)
    result = KC.FourierDiagrams{C5WeakCoeff,Boson,3,1}()
    outer_terms = c5_weak_cut_terms(Γ2, c5_weak_ψ)
    bubble_terms = c5_weak_cut_terms(Γ2, c5_weak_χ)
    component = c5_weak_component(response_wigner, component_kind)

    for outer in outer_terms
        outer.kind === component_kind || continue
        carrier_graph, carrier_kinematic = c5_weak_outer_carrier_graph(outer)
        wigner_graph = c5_weak_wigner_graph(component, carrier_graph, carrier_kinematic)
        KC.response_momentum(wigner_graph) == KC.response_momentum(carrier_graph) ||
            error("C5 changed the routed pair-response momentum")

        atomic = only(c for c in outer.internal if KC.field_family(c.out) == c5_weak_ψ)
        hs = only(c for c in outer.internal if KC.field_family(c.out) == c5_weak_χ)
        response_kind = c5_weak_response_kind(KC.response_component(wigner_graph))
        response_kind === c5_weak_line_kind(hs) ||
            error("C5 changed response Keldysh identity")

        for (Ωsymbol, response_coefficient) in c5_weak_response_terms(response_kind, sector)
            Ωkind = c5_weak_Ωkind(Ωsymbol)
            for bubble in bubble_terms
                bubble.kind === Ωkind || continue
                contractions = KC.Contraction{Boson}[atomic]
                append!(contractions, c5_weak_attach_bubble(hs, bubble))
                KC._twopi_internal_is_physical(contractions) || continue
                sunset, kinematic = c5_weak_anchored_sunset(outer, contractions)
                coefficient = -outer.coefficient * response_coefficient * bubble.coefficient
                KC._push_fourier!(result, sunset, coefficient, kinematic)
            end
        end
    end
    return result
end

function c5_weak_expand_self_energy(Γ2, response_wigner, sector::Symbol, parameter)
    return KC.FourierSelfEnergy{C5WeakCoeff,Boson,2,3,1}(
        c5_weak_expand_component(Γ2, response_wigner, :K, sector),
        c5_weak_expand_component(Γ2, response_wigner, :R, sector),
        c5_weak_expand_component(Γ2, response_wigner, :A, sector),
        parameter,
        c5_weak_ψ,
    )
end

function c5_weak_canonical_collision(ΣF)
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    kinetic = kinetic_expression(ΣW)
    off_shell = off_shell_collision_expression(kinetic)
    return KC.canonical_frequency_collision(spectral_dispersive_collision(off_shell))
end

@testset "C5 weak response collapses to ordinary KC kinetics" begin
    Γ2 = TwoPIEffectiveAction(c5_weak_hs_interaction(), Val(2), Val(3))
    composite = KC.twopi_composite_fourier_self_energy(
        Γ2, c5_weak_ψ, c5_weak_χ; coherent_parameter=:g, loss_parameter=:γ
    )
    response_wigner = KC.response_aware_wigner_transform(composite)
    response_kinetic = kinetic_expression(response_wigner)
    response_collision = off_shell_collision_expression(response_kinetic)

    @test KC.response_polarization(response_kinetic) ===
        KC.response_polarization(response_collision)
    @test KC.target_family(response_collision) === c5_weak_ψ

    sectors = (
        (:g2, KC.ParameterMonomial(:g)^2),
        (:gγ, KC.ParameterMonomial(:g) * KC.ParameterMonomial(:γ)),
        (:γ2, KC.ParameterMonomial(:γ)^2),
    )
    for (sector, parameter) in sectors
        expanded = c5_weak_expand_self_energy(Γ2, response_wigner, sector, parameter)
        direct = c5_weak_direct_skeleton(parameter)
        if sector === :g2
            @test isequal(expanded.keldysh, direct.keldysh)
            @test isequal(expanded.retarded, direct.retarded)
            @test isequal(expanded.advanced, direct.advanced)
        else
            expanded_collision = c5_weak_canonical_collision(expanded)
            direct_collision = c5_weak_canonical_collision(direct)
            @test KC.canonical_frequency_expressions(expanded_collision) ==
                KC.canonical_frequency_expressions(direct_collision)
            if sector === :γ2
                expanded_reduction = KC.reduce_canonical_trotter_frequencies(
                    expanded_collision
                )
                direct_reduction = KC.reduce_canonical_trotter_frequencies(direct_collision)
                @test KC.canonical_trotter_constants(expanded_reduction) ==
                    KC.canonical_trotter_constants(direct_reduction)
                @test KC.canonical_trotter_boundary_expressions(expanded_reduction) ==
                    KC.canonical_trotter_boundary_expressions(direct_reduction)
            end
        end
    end
end
