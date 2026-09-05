##########################################
#       dressed green's function
##########################################
"""
$(DocStringExtensions.TYPEDEF)

A structure representing a dressed propagator in the Retarded-Advanced-Keldysh basis.
Its coefficient representation, statistics, perturbation order, and diagram shape are
encoded in the type.

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
    "Parameters of the perturbation series"
    parameter::CSym
end

function DressedPropagator(
    keldysh::Diagrams{C,S,E1,E2},
    retarded::Diagrams{C,S,E1,E2},
    advanced::Diagrams{C,S,E1,E2},
    ::Val{O},
    parameter::CSym,
) where {C<:Number,S<:Statistics,O,E1,E2}
    return DressedPropagator{C,S,O,E1,E2}(keldysh, retarded, advanced, parameter)
end

order(::DressedPropagator{C,S,O}) where {C,S,O} = O
statistics(::DressedPropagator{C,S}) where {C,S} = S
parameters(d::DressedPropagator) = d.parameter

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

"""
    DressedPropagator(L::InteractionLagrangian, ::Val{order}, ::Val{edges}; target, kwargs...)

For a single field family, the target propagator is inferred. Multi-family interactions
require `target` to select the physical field family.

All the same-coordinate advanced propagators are converted to retarded propagators when
`simplify=true`.
"""
function DressedPropagator(
    L::InteractionLagrangian{C,Boson},
    ::Val{O},
    ::Val{E};
    target=nothing,
    simplify=true,
    _set_reg_to_zero=true,
    kwargs...,
) where {C<:Number,O,E}
    @assert number_of_propagators(L) * O + 1 == E "The supplied Val{edges} must equal the interaction's propagator count times Val{order}, plus the external propagator"
    qfield, cfield = propagator_fields(L, target)

    keldysh = wick_contraction(
        cfield(Out()) * bar(cfield)(In()),
        L,
        Val(O),
        Val(E);
        simplify,
        _set_reg_to_zero,
        kwargs...,
    )
    retarded = wick_contraction(
        cfield(Out()) * bar(qfield)(In()),
        L,
        Val(O),
        Val(E);
        simplify,
        _set_reg_to_zero,
        kwargs...,
    )
    advanced = wick_contraction(
        qfield(Out()) * bar(cfield)(In()),
        L,
        Val(O),
        Val(E);
        simplify,
        _set_reg_to_zero,
        kwargs...,
    )

    for component in (keldysh, retarded, advanced)
        filter_nonzero!(component)
    end

    return DressedPropagator(keldysh, retarded, advanced, Val(O), parameters(L)^O)
end

"""Collection of dressed propagators with distinct perturbation-parameter monomials."""
struct DressedPropagatorSum{K,GS,O}
    arguments::Dict{K,GS}
end

SymbolicUtils.arguments(d::DressedPropagatorSum) = d.arguments
order(::DressedPropagatorSum{K,GS,O}) where {K,GS,O} = O
parameters(d::DressedPropagatorSum) = map(G -> G.parameter, values(arguments(d)))

"""
    DressedPropagator(Ls::LagrangianSum, ::Val{order}, ::Val{edges}; target, kwargs...)
"""
function DressedPropagator(
    Ls::LagrangianSum{C,Boson},
    ::Val{O},
    ::Val{E};
    target=nothing,
    simplify=true,
    _set_reg_to_zero=true,
    kwargs...,
) where {C<:Number,O,E}
    @assert all(number_of_propagators(L) * O + 1 == E for L in arguments(Ls)) "All LagrangianSum terms must produce the supplied number of propagator edges"
    qfield, cfield = propagator_fields(first(arguments(Ls)), target)

    simplify = isa(simplify, Bool) ? fill(simplify, length(Ls)) : simplify

    keldysh_pairs = wick_contraction(
        cfield(Out()) * bar(cfield)(In()),
        Ls,
        Val(O),
        Val(E);
        simplify,
        _set_reg_to_zero,
        kwargs...,
    )
    retarded_pairs = wick_contraction(
        cfield(Out()) * bar(qfield)(In()),
        Ls,
        Val(O),
        Val(E);
        simplify,
        _set_reg_to_zero,
        kwargs...,
    )
    advanced_pairs = wick_contraction(
        qfield(Out()) * bar(cfield)(In()),
        Ls,
        Val(O),
        Val(E);
        simplify,
        _set_reg_to_zero,
        kwargs...,
    )

    D = diagram_coefficient_type(C)
    GS = DressedPropagator{D,Boson,O,E,max_edges(O)}
    dict = Dict{CSym,GS}()
    for idx in eachindex(keldysh_pairs)
        components = last.((keldysh_pairs[idx], retarded_pairs[idx], advanced_pairs[idx]))
        for component in components
            filter_nonzero!(component)
            _simplify_prefactors!(component)
        end
        parameter = first(keldysh_pairs[idx])
        dict[parameter] = DressedPropagator(components..., Val(O), parameter)
    end

    return DressedPropagatorSum{CSym,GS,O}(dict)
end
