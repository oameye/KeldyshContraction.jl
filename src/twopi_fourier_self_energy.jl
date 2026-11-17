"""Return the physical external propagator types for a bosonic self-energy component."""
function _twopi_fourier_external_types(component::UInt8)
    if component == 0x01
        return PropagatorType.Retarded, PropagatorType.Keldysh
    elseif component == 0x02
        return PropagatorType.Keldysh, PropagatorType.Advanced
    elseif component == 0x03
        return PropagatorType.Retarded, PropagatorType.Advanced
    end
    return throw(
        ArgumentError("structural cl-cl self-energy has no physical external legs")
    )
end

function _twopi_fourier_out_anchor(
    cut::Contraction{Boson}, target::FieldFamily{Boson}, desired::PropagatorType.T
)
    for index in (Classical, Quantum)
        candidate = Contraction(target[index](Out()), cut.in)
        is_qq_contraction(candidate) && continue
        propagator_type(candidate) === desired && return candidate
    end
    return error("could not construct the requested outgoing 2PI Fourier anchor")
end

function _twopi_fourier_in_anchor(
    cut::Contraction{Boson}, target::FieldFamily{Boson}, desired::PropagatorType.T
)
    for index in (Classical, Quantum)
        candidate = Contraction(cut.out, bar(target[index])(In()))
        is_qq_contraction(candidate) && continue
        propagator_type(candidate) === desired && return candidate
    end
    return error("could not construct the requested incoming 2PI Fourier anchor")
end

function _twopi_fourier_anchored_diagram(
    contractions::FixedVector{E,Contraction{Boson}},
    cut_index::Int,
    target::FieldFamily{Boson},
    component::UInt8,
    ::Val{ST},
) where {E,ST}
    cut = contractions[cut_index]
    out_type, in_type = _twopi_fourier_external_types(component)

    anchored = Contraction{Boson}[]
    sizehint!(anchored, E + 1)
    for index in eachindex(contractions)
        index == cut_index && continue
        push!(anchored, contractions[index])
    end
    push!(anchored, _twopi_fourier_out_anchor(cut, target, out_type))
    push!(anchored, _twopi_fourier_in_anchor(cut, target, in_type))
    return Diagram(anchored, Val(E + 1), Val(ST))
end

"""
    twopi_fourier_self_energy(Γ₂, target)

Lower a bosonic [`TwoPIEffectiveAction`](@ref) directly to the exact Fourier-space
`FourierSelfEnergy` representation for `target`.

The target full-propagator line is cut while its two interaction-vertex endpoints are still known.
Physical external propagator segments are attached to those endpoints, the complete two-point graph
is Fourier-routed and derivative-lowered with KC's ordinary Fourier machinery, and only then are the
external segments amputated. This ordering preserves external momentum factors and is therefore
intentionally different from `fourier_transform(SelfEnergy(Γ₂, target))`, which is unsafe and
remains rejected.

The returned `FourierSelfEnergy` uses the ordinary downstream representation and can be passed
directly to `wigner_transform`. The operation is purely symbolic; it does not solve the Dyson or
Kadanoff--Baym equations.
"""
function twopi_fourier_self_energy(
    Γ::TwoPIEffectiveAction{C,Boson,O,E,E2}, target::FieldFamily{Boson}
) where {C<:Number,O,E,E2}
    target in field_families(Γ) ||
        throw(ArgumentError("target field family is absent from the 2PI effective action"))

    SE = E - 1
    ST = max_edges(O)
    D = FourierDiagrams{C,Boson,SE,ST}
    self_energy = SmallCollections.SmallDict{3,PropagatorType.T,D}((
        PropagatorType.Advanced => D(),
        PropagatorType.Retarded => D(),
        PropagatorType.Keldysh => D(),
    ))
    derivative_factor = _twopi_derivative_factor(Boson, C)

    for (vacuum, coefficient) in Γ
        contractions = twopi_contractions(vacuum)
        for cut_index in eachindex(contractions)
            cut = contractions[cut_index]
            isequal(field_family(cut.out), target) || continue

            component = _twopi_physical_component(
                keldysh_index(cut.in), keldysh_index(cut.out)
            )
            iszero(component) && continue

            internal = Contraction{Boson}[
                contractions[index] for
                index in eachindex(contractions) if index != cut_index
            ]
            _twopi_internal_is_physical(internal) || continue

            source = _twopi_fourier_anchored_diagram(
                contractions, cut_index, target, component, Val(ST)
            )
            routed, kinematic, _ = _canonical_fourier_source(source)
            amputated = _amputate_fourier_graph(routed, Val(SE), Val(ST))

            destination = if component == 0x01
                self_energy[PropagatorType.Retarded]
            elseif component == 0x02
                self_energy[PropagatorType.Advanced]
            else
                self_energy[PropagatorType.Keldysh]
            end
            _push_fourier!(
                destination, amputated, derivative_factor * coefficient, kinematic
            )
        end
    end

    return FourierSelfEnergy{C,Boson,O,SE,ST}(
        self_energy[PropagatorType.Keldysh],
        self_energy[PropagatorType.Retarded],
        self_energy[PropagatorType.Advanced],
        parameters(Γ),
        target,
    )
end
