# Direct self-energy construction. This bypasses materialization of reducible two-point diagrams
# while preserving the existing `SelfEnergy(DressedPropagator(...))` path as an independent
# compatibility oracle.

function _direct_self_energy_from_keldysh(
    keldysh::Diagrams{C,Boson,E,ST},
    ::Val{O},
    parameter::ParameterMonomial,
    target::FieldFamily{Boson},
) where {C<:Number,E,ST,O}
    SE = E - 2
    D = Diagrams{C,Boson,SE,ST}
    self_energy = SmallCollections.SmallDict{3,PropagatorType.T,D}((
        PropagatorType.Advanced => D(),
        PropagatorType.Retarded => D(),
        PropagatorType.Keldysh => D(),
    ))

    # This deliberately mirrors `_self_energy(G)`: all three self-energy components are classified
    # from the physical contractions in the Keldysh two-point function. The separate retarded and
    # advanced dressed-propagator components are not inputs to self-energy extraction.
    construct_self_energy!(self_energy, keldysh)

    _simplify_prefactors!(self_energy[PropagatorType.Keldysh])
    _simplify_prefactors!(self_energy[PropagatorType.Retarded])
    _simplify_prefactors!(self_energy[PropagatorType.Advanced])

    return SelfEnergy{C,Boson,O,SE,ST}(
        self_energy[PropagatorType.Keldysh],
        self_energy[PropagatorType.Retarded],
        self_energy[PropagatorType.Advanced],
        parameter,
        target,
    )
end

function _direct_self_energy(
    L::InteractionLagrangian{C,Boson},
    order::Val{O},
    edges::Val{E};
    target=nothing,
    simplify=true,
    preserve_regularisation=false,
) where {C<:Number,O,E}
    @assert number_of_propagators(L) * O + 1 == E "The supplied Val{edges} must equal the interaction's propagator count times Val{order}, plus the external propagator"

    fields = propagator_fields(L, target)
    selected_target = field_family(first(fields))
    keldysh_external = first(propagator_external_products(Boson, fields...))
    set_reg_to_zero = !preserve_regularisation

    keldysh = _onepi_wick_contraction(
        keldysh_external, L, order, edges; simplify, _set_reg_to_zero=set_reg_to_zero
    )
    filter_nonzero!(keldysh)

    return _direct_self_energy_from_keldysh(
        keldysh, order, parameters(L)^O, selected_target
    )
end

"""
    SelfEnergy(L::InteractionLagrangian, ::Val{order}, ::Val{edges};
        target=nothing, simplify=true, preserve_regularisation=false)

Construct the fixed-order bosonic 1PI self-energy directly from an interaction without first
materializing the reducible dressed two-point diagrams. The result is semantically identical to
`SelfEnergy(DressedPropagator(L, order, edges; ...))`; the latter remains available when the full
dressed propagator is needed independently.
"""
function SelfEnergy(
    L::InteractionLagrangian{C,Boson}, order::Val{O}, edges::Val{E}; kwargs...
) where {C<:Number,O,E}
    return _direct_self_energy(L, order, edges; kwargs...)
end
