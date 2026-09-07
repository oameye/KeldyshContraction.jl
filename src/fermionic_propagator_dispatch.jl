# Fermionic high-level entry points stay statistics-specific so inference preserves the
# relation between the external Grassmann products and the interaction statistics. The
# diagram-generation implementations themselves are shared with the bosonic path.
function DressedPropagator(
    L::InteractionLagrangian{C,Fermion},
    order::Val{O},
    edges::Val{E};
    target=nothing,
    simplify=true,
    _set_reg_to_zero=true,
    kwargs...,
) where {C<:Number,O,E}
    fields = propagator_fields(L, target)
    external_products = propagator_external_products(Fermion, fields...)
    return _dressed_propagator(
        L, external_products, order, edges; simplify, _set_reg_to_zero, kwargs...
    )
end

function DressedPropagator(
    Ls::LagrangianSum{C,Fermion},
    order::Val{O},
    edges::Val{E};
    target=nothing,
    simplify=true,
    _set_reg_to_zero=true,
    kwargs...,
) where {C<:Number,O,E}
    fields = propagator_fields(first(arguments(Ls)), target)
    external_products = propagator_external_products(Fermion, fields...)
    return _dressed_propagator_sum(
        Ls, external_products, order, edges; simplify, _set_reg_to_zero, kwargs...
    )
end
