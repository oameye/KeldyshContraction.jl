##########################################
#       dressed green's function
##########################################
"""
$(DocStringExtensions.TYPEDEF)

A structure representing a dressed propagator in the Retarded-Advanced-Keldysh basis.
Its coefficient representation, statistics, perturbation order, and diagram shape are
encoded in the type. The selected external physical field family is retained explicitly.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct DressedPropagator{C<:Number,S<:Statistics,O,E1,E2}
    "The Keldysh component of the propagator"
    keldysh::Diagrams{C,S,E1,E2}
    "The retarded component of the propagator"
    retarded::Diagrams{C,S,E1,E2}
    "The advanced component of the propagator"
    advanced::Diagrams{C,S,E1,E2}
    "Canonical perturbation-parameter monomial"
    parameter::ParameterMonomial
    "Physical field family selected for the external two-point function"
    target::FieldFamily{S}
end

function DressedPropagator(
    keldysh::Diagrams{C,S,E1,E2},
    retarded::Diagrams{C,S,E1,E2},
    advanced::Diagrams{C,S,E1,E2},
    ::Val{O},
    parameter::ParameterMonomial,
    target::FieldFamily{S},
) where {C<:Number,S<:Statistics,O,E1,E2}
    return DressedPropagator{C,S,O,E1,E2}(keldysh, retarded, advanced, parameter, target)
end

function Base.isequal(
    a::DressedPropagator{C,S,O,E1,E2}, b::DressedPropagator{C,S,O,E1,E2}
) where {C,S,O,E1,E2}
    return isequal(a.keldysh, b.keldysh) &&
           isequal(a.retarded, b.retarded) &&
           isequal(a.advanced, b.advanced) &&
           isequal(a.parameter, b.parameter) &&
           isequal(a.target, b.target)
end
function Base.:(==)(
    a::DressedPropagator{C,S,O,E1,E2}, b::DressedPropagator{C,S,O,E1,E2}
) where {C,S,O,E1,E2}
    return isequal(a, b)
end
function Base.hash(d::DressedPropagator, h::UInt)
    return hash((d.keldysh, d.retarded, d.advanced, d.parameter, d.target), h)
end

order(::DressedPropagator{C,S,O}) where {C,S,O} = O
statistics(::DressedPropagator{C,S}) where {C,S} = S
parameters(d::DressedPropagator) = d.parameter
target_family(d::DressedPropagator) = d.target

function structural_zero(
    ::Type{Boson}, ::Type{Diagrams{C,Boson,E1,E2}}
) where {C<:Number,E1,E2}
    return Diagrams{C,Boson,E1,E2}()
end

"""
    matrix(G::DressedPropagator)

Return the bosonic Retarded-Advanced-Keldysh matrix
```math
\\hat{G}=\\begin{pmatrix}G^K&G^R\\\\G^A&0\\end{pmatrix}.
```
"""
function matrix(G::DressedPropagator{C,Boson,O,E1,E2}) where {C<:Number,O,E1,E2}
    return matrix(Boson, G)
end
function matrix(
    ::Type{Boson}, G::DressedPropagator{C,Boson,O,E1,E2}
) where {C<:Number,O,E1,E2}
    D = Diagrams{C,Boson,E1,E2}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = G.keldysh
    result[1, 2] = G.retarded
    result[2, 1] = G.advanced
    result[2, 2] = structural_zero(Boson, D)
    return result
end

function propagator_fields(L::InteractionLagrangian{C,Boson}, ::Nothing) where {C<:Number}
    family = target_family(L)
    return family[Quantum], family[Classical]
end
function propagator_fields(
    L::InteractionLagrangian{C,Boson}, target::FieldFamily{Boson}
) where {C<:Number}
    family = target_family(L, target)
    return family[Quantum], family[Classical]
end

function propagator_external_products(
    ::Type{Boson}, qfield::Field{Boson}, cfield::Field{Boson}
)
    return (
        cfield(Out()) * bar(cfield)(In()),
        cfield(Out()) * bar(qfield)(In()),
        qfield(Out()) * bar(cfield)(In()),
    )
end

function _dressed_propagator(
    L::InteractionLagrangian{C,S},
    external_products::Tuple{QMul{P1,S},QMul{P2,S},QMul{P3,S}},
    target::FieldFamily{S},
    ::Val{O},
    ::Val{E};
    simplify=true,
    _set_reg_to_zero=true,
    kwargs...,
) where {C<:Number,S<:Statistics,P1<:Number,P2<:Number,P3<:Number,O,E}
    @assert number_of_propagators(L) * O + 1 == E "The supplied Val{edges} must equal the interaction's propagator count times Val{order}, plus the external propagator"
    keldysh_inout, retarded_inout, advanced_inout = external_products

    keldysh = wick_contraction(
        keldysh_inout, L, Val(O), Val(E); simplify, _set_reg_to_zero, kwargs...
    )
    retarded = wick_contraction(
        retarded_inout, L, Val(O), Val(E); simplify, _set_reg_to_zero, kwargs...
    )
    advanced = wick_contraction(
        advanced_inout, L, Val(O), Val(E); simplify, _set_reg_to_zero, kwargs...
    )

    for component in (keldysh, retarded, advanced)
        filter_nonzero!(component)
    end

    return DressedPropagator(keldysh, retarded, advanced, Val(O), parameters(L)^O, target)
end

"""
    DressedPropagator(L::InteractionLagrangian, ::Val{order}, ::Val{edges};
        target=nothing, simplify=true, preserve_regularisation=false, kwargs...)

For a single field family, the target propagator is inferred. Multi-family interactions
require `target` to select the physical field family. The selected family is retained in the
returned propagator. Statistics dispatch selects the external K/R/A field products; diagram
generation then shares one implementation.

All the same-coordinate advanced propagators are converted to retarded propagators when
`simplify=true`. Set `preserve_regularisation=true` when equal-time/Trotter shifts must
remain attached to the generated diagrams; the default removes those shifts after the Wick
contraction, preserving the historical high-level behavior.
"""
function DressedPropagator(
    L::InteractionLagrangian{C,Boson},
    order::Val{O},
    edges::Val{E};
    target=nothing,
    simplify=true,
    preserve_regularisation=false,
    kwargs...,
) where {C<:Number,O,E}
    fields = propagator_fields(L, target)
    selected_target = field_family(first(fields))
    external_products = propagator_external_products(Boson, fields...)
    return _dressed_propagator(
        L,
        external_products,
        selected_target,
        order,
        edges;
        simplify,
        _set_reg_to_zero=!preserve_regularisation,
        kwargs...,
    )
end

"""Collection of dressed propagators with distinct perturbation-parameter monomials."""
struct DressedPropagatorSum{GS,O}
    arguments::Dict{ParameterMonomial,GS}
end

SymbolicUtils.arguments(d::DressedPropagatorSum) = d.arguments
function Base.getindex(d::DressedPropagatorSum, parameter)
    return d.arguments[parameter_monomial(parameter)]
end
order(::DressedPropagatorSum{GS,O}) where {GS,O} = O
parameters(d::DressedPropagatorSum) = collect(keys(d.arguments))

function _dressed_propagator_sum(
    Ls::LagrangianSum{C,S},
    external_products::Tuple{QMul{P1,S},QMul{P2,S},QMul{P3,S}},
    target::FieldFamily{S},
    ::Val{O},
    ::Val{E};
    simplify=true,
    _set_reg_to_zero=true,
    kwargs...,
) where {C<:Number,S<:Statistics,P1<:Number,P2<:Number,P3<:Number,O,E}
    @assert all(number_of_propagators(L) * O + 1 == E for L in arguments(Ls)) "All LagrangianSum terms must produce the supplied number of propagator edges"
    keldysh_inout, retarded_inout, advanced_inout = external_products
    simplify_flags = isa(simplify, Bool) ? fill(simplify, length(Ls)) : simplify

    keldysh_pairs = wick_contraction(
        keldysh_inout,
        Ls,
        Val(O),
        Val(E);
        simplify=simplify_flags,
        _set_reg_to_zero,
        kwargs...,
    )
    retarded_pairs = wick_contraction(
        retarded_inout,
        Ls,
        Val(O),
        Val(E);
        simplify=simplify_flags,
        _set_reg_to_zero,
        kwargs...,
    )
    advanced_pairs = wick_contraction(
        advanced_inout,
        Ls,
        Val(O),
        Val(E);
        simplify=simplify_flags,
        _set_reg_to_zero,
        kwargs...,
    )

    D = diagram_coefficient_type(C)
    GS = DressedPropagator{D,S,O,E,max_edges(O)}
    dict = Dict{ParameterMonomial,GS}()
    for idx in eachindex(keldysh_pairs)
        components = last.((keldysh_pairs[idx], retarded_pairs[idx], advanced_pairs[idx]))
        for component in components
            filter_nonzero!(component)
            _simplify_prefactors!(component)
        end
        parameter = first(keldysh_pairs[idx])
        dict[parameter] = DressedPropagator(components..., Val(O), parameter, target)
    end

    return DressedPropagatorSum{GS,O}(dict)
end

"""
    DressedPropagator(Ls::LagrangianSum, ::Val{order}, ::Val{edges};
        target=nothing, simplify=true, preserve_regularisation=false, kwargs...)

Construct a perturbative propagator for a sum of interactions with common field families.
Statistics dispatch selects the R/A/K external products while diagram generation and
parameter-monomial accumulation share one implementation. The selected external family is
retained by every propagator in the sum. Set `preserve_regularisation=true` to retain
finite equal-time/Trotter shifts in every generated component.
"""
function DressedPropagator(
    Ls::LagrangianSum{C,Boson},
    order::Val{O},
    edges::Val{E};
    target=nothing,
    simplify=true,
    preserve_regularisation=false,
    kwargs...,
) where {C<:Number,O,E}
    fields = propagator_fields(first(arguments(Ls)), target)
    selected_target = field_family(first(fields))
    external_products = propagator_external_products(Boson, fields...)
    return _dressed_propagator_sum(
        Ls,
        external_products,
        selected_target,
        order,
        edges;
        simplify,
        _set_reg_to_zero=!preserve_regularisation,
        kwargs...,
    )
end
