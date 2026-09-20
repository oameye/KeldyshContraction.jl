# Private direct self-energy construction. This bypasses materialization of reducible two-point
# diagrams while preserving the existing public `SelfEnergy(DressedPropagator(...))` path as an
# independent compatibility oracle.

function _direct_self_energy_from_components(
    keldysh::Diagrams{C,Boson,E,ST},
    retarded::Diagrams{C,Boson,E,ST},
    advanced::Diagrams{C,Boson,E,ST},
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

    construct_self_energy!(self_energy, keldysh)
    construct_self_energy!(self_energy, retarded)
    construct_self_energy!(self_energy, advanced)

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
    external_products = propagator_external_products(Boson, fields...)
    set_reg_to_zero = !preserve_regularisation

    components = map(external_products) do external_product
        diagrams = _onepi_wick_contraction(
            external_product, L, order, edges; simplify, _set_reg_to_zero=set_reg_to_zero
        )
        filter_nonzero!(diagrams)
        return diagrams
    end

    return _direct_self_energy_from_components(
        components..., order, parameters(L)^O, selected_target
    )
end
