using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields weak_coupling_ψ::Boson weak_coupling_χ::Boson

const WeakCouplingCoeff = Complex{Rational{Int}}

function weak_coupling_microscopic_interaction()
    c = weak_coupling_ψ[Classical]
    q = weak_coupling_ψ[Quantum]
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

function weak_coupling_hs_interaction()
    ψc = weak_coupling_ψ[Classical]
    ψq = weak_coupling_ψ[Quantum]
    χc = weak_coupling_χ[Classical]
    χq = weak_coupling_χ[Quantum]

    # KC uses the symmetric Keldysh rotation ψ±=(ψc±ψq)/√2 (and likewise for χ).
    # Rotating the contour HS vertex
    #   -(i/√2) [ψ₊² χ̄₊ - ψ₋² χ̄₋ + h.c.]
    # directly into this convention gives the physical cubic vertex below. No extra
    # h=1/√2 normalization remains after the rotation.
    forward = (1 // 2) * ψc^2 * bar(χq) + ψc * ψq * bar(χc) + (1 // 2) * ψq^2 * bar(χq)
    return ChargedInteractionLagrangian(
        -im * (forward + bar(forward)),
        weak_coupling_ψ => 1,
        weak_coupling_χ => 2;
        parameter=:h,
    )
end

function weak_coupling_microscopic_coordinate_skeleton(parameter)
    L = weak_coupling_microscopic_interaction()
    G = DressedPropagator(L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false)
    return KC.skeleton_self_energy(G[parameter])
end

function weak_coupling_microscopic_skeleton(parameter)
    L = weak_coupling_microscopic_interaction()
    G = DressedPropagator(L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false)
    GF = fourier_transform(G[parameter])
    return KC.skeleton_self_energy(GF)
end

function weak_coupling_microscopic_full_self_energy(parameter)
    L = weak_coupling_microscopic_interaction()
    G = DressedPropagator(L, Val(2), Val(5); simplify=true, _set_reg_to_zero=false)
    GF = fourier_transform(G[parameter])
    return SelfEnergy(GF)
end

function weak_coupling_line_kind(edge)
    KC.is_keldysh(edge) && return :K
    KC.is_retarded(edge) && return :R
    KC.is_advanced(edge) && return :A
    return error("unexpected bosonic propagator type")
end

function weak_coupling_component(Σ, kind::Symbol)
    kind === :K && return KC.keldysh_component(Σ)
    kind === :R && return KC.retarded_component(Σ)
    kind === :A && return KC.advanced_component(Σ)
    return error("unexpected Keldysh component")
end

function weak_coupling_component_kind(component::UInt8)
    component == 0x01 && return :R
    component == 0x02 && return :A
    component == 0x03 && return :K
    return error("unexpected physical 2PI component")
end

function weak_coupling_component_census(component)
    records = NamedTuple[]
    for (graph, contributions) in component
        signature = sort!(
            collect(
                weak_coupling_line_kind(edge) for edge in KC.contractions(graph.coordinate)
            ),
        )
        for contribution in contributions
            push!(
                records,
                (;
                    signature=Tuple(signature),
                    coefficient=contribution.coefficient,
                    routing=KC.edge_momenta(graph),
                ),
            )
        end
    end
    sort!(records; by=repr)
    return records
end

function weak_coupling_response_terms(kind::Symbol, sector::Symbol)
    one_c = one(WeakCouplingCoeff)
    i_c = complex(0 // 1, 1 // 1)

    # In KC's symmetric cl/q normalization the bare HS matrix is
    #   D₀ = [-2γ  -iλ; -iλ̄  0].
    # Therefore Dreg^(2)=D₀ Ω D₀ keeps the C2 R/A coefficients while the terms
    # involving the Keldysh bare contact are doubled relative to the asymmetric
    # source convention used in the original note.
    if kind === :R
        coefficient = if sector === :g2
            -one_c
        elseif sector === :gγ
            2 * i_c
        else
            one_c
        end
        return ((:R, coefficient),)
    elseif kind === :A
        coefficient = if sector === :g2
            -one_c
        elseif sector === :gγ
            -2 * i_c
        else
            one_c
        end
        return ((:A, coefficient),)
    elseif kind === :K
        sector === :g2 && return ((:K, -one_c),)
        sector === :gγ && return ((:R, 2 * i_c), (:A, 2 * i_c))
        sector === :γ2 && return ((:K, -one_c), (:R, 2 * one_c), (:A, -2 * one_c))
    end
    return error("unexpected response sector")
end

# Keep the original cut line together with its two interaction vertices. This is the
# provenance that an ordinary amputated SelfEnergy intentionally no longer stores.
function weak_coupling_twopi_cut_terms(Γ2, target)
    records = NamedTuple[]
    derivative_factor = convert(WeakCouplingCoeff, im)
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
                    kind=weak_coupling_component_kind(component),
                    coefficient=derivative_factor *
                                convert(WeakCouplingCoeff, vacuum_coefficient),
                ),
            )
        end
    end
    return records
end

function weak_coupling_reposition(field, row_from, row_to, column_from, column_to)
    p = KC.position(field)
    p == row_from && return field(row_to)
    p == column_from && return field(column_to)
    return error("pair bubble contains a vertex unrelated to its χ cut")
end

function weak_coupling_attach_bubble(outer_hs, bubble)
    # The χ derivative is transposed: cut.in is the response row and cut.out its column.
    # D₀ Ω D₀ therefore identifies cut.in with D.out and cut.out with D.in. Keep the
    # compiler-generated pair lines themselves untouched; endpoint attachment supplies
    # the required coordinate ordering, while their Keldysh labels remain those of Ω.
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
                weak_coupling_reposition(
                    contraction.out, row_from, row_to, column_from, column_to
                ),
                weak_coupling_reposition(
                    contraction.in, row_from, row_to, column_from, column_to
                ),
            ),
        )
    end
    return attached
end

function weak_coupling_anchored_fourier_sunset(outer, contractions)
    # Reuse the C0 ordering exactly: physical external ψ segments stay attached during
    # canonical Fourier routing and are amputated only afterwards. Marker-only anchoring
    # followed by a second `fourier_transform` is insufficient because that transform
    # recanonicalizes the stripped sunset and can exchange the two self-energy vertices.
    length(contractions) == 3 ||
        error("anchored response substitution produced an invalid sunset")
    component = KC._twopi_physical_component(
        KC.keldysh_index(outer.cut.in), KC.keldysh_index(outer.cut.out)
    )
    out_type, in_type = KC._twopi_fourier_external_types(component)
    out_anchor = KC._twopi_fourier_out_anchor(outer.cut, weak_coupling_ψ, out_type)
    in_anchor = KC._twopi_fourier_in_anchor(outer.cut, weak_coupling_ψ, in_type)
    source = KC.Diagram(
        KC.Contraction{Boson}[contractions..., out_anchor, in_anchor], Val(5), Val(1)
    )
    routed, kinematic, _ = KC._canonical_fourier_source(source)
    amputated = KC._amputate_fourier_graph(routed, Val(3), Val(1))
    return amputated, kinematic
end

function weak_coupling_hs_expanded_component(Γ2, component_kind::Symbol, sector::Symbol)
    result = KC.FourierDiagrams{WeakCouplingCoeff,Boson,3,1}()
    outer_terms = weak_coupling_twopi_cut_terms(Γ2, weak_coupling_ψ)
    bubble_terms = weak_coupling_twopi_cut_terms(Γ2, weak_coupling_χ)

    for outer in outer_terms
        outer.kind === component_kind || continue
        atomic = only(
            contraction for contraction in outer.internal if
            KC.field_family(contraction.out) == weak_coupling_ψ
        )
        hs = only(
            contraction for contraction in outer.internal if
            KC.field_family(contraction.out) == weak_coupling_χ
        )
        hs_kind = weak_coupling_line_kind(hs)

        for (Ωkind, response_coefficient) in weak_coupling_response_terms(hs_kind, sector)
            for bubble in bubble_terms
                bubble.kind === Ωkind || continue
                contractions = KC.Contraction{Boson}[atomic]
                append!(contractions, weak_coupling_attach_bubble(hs, bubble))
                KC._twopi_internal_is_physical(contractions) || continue
                sunset, kinematic = weak_coupling_anchored_fourier_sunset(
                    outer, contractions
                )

                # Each compiler-generated O(h²) 2PI cut carries KC's common -i kernel phase.
                # A physical response insertion composes two such kernels, so relative to the
                # directly generated O(g²,gγ,γ²) quartic sunset the exact Wick phase is i²=-1.
                coefficient = -outer.coefficient * response_coefficient * bubble.coefficient
                KC._push_fourier!(result, sunset, coefficient, kinematic)
            end
        end
    end
    return result
end

function weak_coupling_hs_fourier_self_energy(Γ2, sector::Symbol, parameter)
    return KC.FourierSelfEnergy{WeakCouplingCoeff,Boson,2,3,1}(
        weak_coupling_hs_expanded_component(Γ2, :K, sector),
        weak_coupling_hs_expanded_component(Γ2, :R, sector),
        weak_coupling_hs_expanded_component(Γ2, :A, sector),
        parameter,
        weak_coupling_ψ,
    )
end

function weak_coupling_canonical_collision(ΣF)
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    KΣ = kinetic_expression(ΣW)
    off_shell = off_shell_collision_expression(KΣ)
    spectral = spectral_dispersive_collision(off_shell)
    return KC.canonical_frequency_collision(spectral)
end

function weak_coupling_statistical_families(collision)
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

function weak_coupling_support_signature(supports)
    return Tuple(
        (
            KC.affine_support_rank(support),
            KC.affine_support_count(support),
            Tuple(Tuple(row) for row in support.dependency_rows),
        ) for support in supports
    )
end

@testset "microscopic second-order skeleton is the C3 Fourier oracle" begin
    sectors = (
        (:g2, KC.ParameterMonomial(:g)^2),
        (:gγ, KC.ParameterMonomial(:g) * KC.ParameterMonomial(:γ)),
        (:γ2, KC.ParameterMonomial(:γ)^2),
    )

    for (label, parameter) in sectors
        Σ = weak_coupling_microscopic_skeleton(parameter)
        @test parameters(Σ) == parameter
        @test target_family(Σ) === weak_coupling_ψ

        for component in (Σ.keldysh, Σ.retarded, Σ.advanced)
            @test !isempty(component)
            for (graph, contributions) in component
                @test graph.external_count == 1
                @test graph.loop_count == 2
                @test collect(KC.topology(graph.coordinate)) == [3]
                @test length(KC.contractions(graph.coordinate)) == 3
                @test all(
                    edge -> KC.field_family(edge.out) == weak_coupling_ψ,
                    KC.contractions(graph.coordinate),
                )
                @test !isempty(contributions)
            end
        end

        @info "C3 microscopic Fourier-sunset census" label keldysh = weak_coupling_component_census(
            Σ.keldysh
        ) retarded = weak_coupling_component_census(Σ.retarded) advanced = weak_coupling_component_census(
            Σ.advanced
        )
    end
end

@testset "KC-native weak-coupling HS response composes into physical sunsets" begin
    Γ2 = TwoPIEffectiveAction(weak_coupling_hs_interaction(), Val(2), Val(3))

    sectors = (
        (:g2, KC.ParameterMonomial(:g)^2),
        (:gγ, KC.ParameterMonomial(:g) * KC.ParameterMonomial(:γ)),
        (:γ2, KC.ParameterMonomial(:γ)^2),
    )

    for (sector, parameter) in sectors
        microscopic = weak_coupling_microscopic_skeleton(parameter)
        for component_kind in (:K, :R, :A)
            hs = weak_coupling_hs_expanded_component(Γ2, component_kind, sector)
            direct = weak_coupling_component(microscopic, component_kind)

            @test !isempty(hs)
            for (graph, contributions) in hs
                @test graph.external_count == 1
                @test graph.loop_count == 2
                @test collect(KC.topology(graph.coordinate)) == [3]
                @test all(
                    edge -> KC.field_family(edge.out) == weak_coupling_ψ,
                    KC.contractions(graph.coordinate),
                )
                @test !isempty(contributions)
            end

            sector === :g2 && @test isequal(hs, direct)
        end
    end

    mixed_parameter = KC.ParameterMonomial(:g) * KC.ParameterMonomial(:γ)
    mixed_hs = weak_coupling_canonical_collision(
        weak_coupling_hs_fourier_self_energy(Γ2, :gγ, mixed_parameter)
    )
    mixed_direct = weak_coupling_canonical_collision(
        weak_coupling_microscopic_skeleton(mixed_parameter)
    )

    @test isempty(KC.shifted_frequency_terms(mixed_hs))
    @test isempty(KC.shifted_frequency_terms(mixed_direct))
    @test KC.canonical_frequency_expressions(mixed_hs) ==
        KC.canonical_frequency_expressions(mixed_direct)
    @test all(==(weak_coupling_ψ), weak_coupling_statistical_families(mixed_hs))

    gamma2_parameter = KC.ParameterMonomial(:γ)^2
    gamma2_hs = weak_coupling_canonical_collision(
        weak_coupling_hs_fourier_self_energy(Γ2, :γ2, gamma2_parameter)
    )
    gamma2_direct = weak_coupling_canonical_collision(
        weak_coupling_microscopic_skeleton(gamma2_parameter)
    )
    gamma2_full = weak_coupling_canonical_collision(
        weak_coupling_microscopic_full_self_energy(gamma2_parameter)
    )

    # The ordinary second-order self-energy retains 32 explicit finite-Trotter terms,
    # reproducing the independent loss oracle. The 2PI skeleton projection removes those
    # non-skeleton contributions before C3 comparison; neither the microscopic skeleton nor
    # the composite HS response therefore carries a shifted statistical degree of freedom.
    @test length(KC.shifted_frequency_terms(gamma2_full)) == 32
    @test isempty(KC.shifted_frequency_terms(gamma2_hs))
    @test isempty(KC.shifted_frequency_terms(gamma2_direct))
    @test all(==(weak_coupling_ψ), weak_coupling_statistical_families(gamma2_hs))

    gamma2_hs_reduction = KC.reduce_canonical_trotter_frequencies(gamma2_hs)
    gamma2_direct_reduction = KC.reduce_canonical_trotter_frequencies(gamma2_direct)

    @test isempty(KC.unresolved_trotter_states(gamma2_hs_reduction))
    @test isempty(KC.unresolved_trotter_states(gamma2_direct_reduction))
    @test KC.canonical_trotter_constants(gamma2_hs_reduction) ==
        KC.canonical_trotter_constants(gamma2_direct_reduction)
    @test KC.canonical_trotter_boundary_expressions(gamma2_hs_reduction) ==
        KC.canonical_trotter_boundary_expressions(gamma2_direct_reduction)

    hs_blockers = KC.blocked_trotter_contributions(gamma2_hs_reduction)
    direct_blockers = KC.blocked_trotter_contributions(gamma2_direct_reduction)
    @test Set(keys(hs_blockers)) == Set(keys(direct_blockers))
    for (key, hs_blocked) in hs_blockers
        direct_blocked = direct_blockers[key]
        @test KC.trotter_frequency_blocker_kind(hs_blocked) ===
            KC.trotter_frequency_blocker_kind(direct_blocked)
        @test KC.trotter_frequency_blocked_index(hs_blocked) ==
            KC.trotter_frequency_blocked_index(direct_blocked)
        @test KC.trotter_frequency_blocked_expression(hs_blocked) ==
            KC.trotter_frequency_blocked_expression(direct_blocked)
        @test weak_coupling_support_signature(
            KC.trotter_frequency_blocked_supports(hs_blocked)
        ) == weak_coupling_support_signature(
            KC.trotter_frequency_blocked_supports(direct_blocked)
        )
    end
end
