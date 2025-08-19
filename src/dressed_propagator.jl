##########################################
#       dressed green's function
##########################################
"""
$(DocStringExtensions.TYPEDEF)

A structure representing dressed propagator in the Retarded-Advanced-Keldysh basis
([`PropagatorType`](@ref)).

# Fields
$(DocStringExtensions.FIELDS)

# Constructor
$(DocStringExtensions.TYPEDSIGNATURES)

Constructs a `DressedPropagator` with the given Keldysh, retarded, and advanced components.
"""
struct DressedPropagator{E1,E2}
    "The Keldysh component of the propagator"
    keldysh::Diagrams{E1,E2}
    "The retarded component of the propagator"
    retarded::Diagrams{E1,E2}
    "The advanced component of the propagator"
    advanced::Diagrams{E1,E2}
    "The order of the dressed propagator in the perturbation series"
    order::Int64
    "Parameters of the perturbation series of the `order`, i.e., ``g^2`` means `order=2`."
    parameter::SymbolicUtils.BasicSymbolic{Number}
end

"""
    matrix(G::DressedPropagator)

Returns the matrix representation of the dressed propagator `G`
in the Retarded-Advanced-Keldysh basis.
```math
\\hat{G}\\left(x_1, x_2\\right)
=\\left(
\\begin{array}{cc}
G^K\\left(x_1, x_2\\right) & G^R\\left(x_1, x_2\\right) \\\\
G^A\\left(x_1, x_2\\right) & 0
\\end{array}
\\right)
```
"""
function matrix(G::DressedPropagator{E1,E2}) where {E1,E2}
    return Diagrams[G.retarded G.keldysh; G.advanced Diagrams{E1,E2}()]
end

"""
    wick_contraction(L::InteractionLagrangian, order=1)


All the same coordinate advanced propagators are converted to retarded propagators.
"""
function DressedPropagator(
    L::InteractionLagrangian,
    order=1;
    simplify=(!should_regularise(L.lagrangian)),
    kwargs...,
)
    ϕ = L.qfield
    ψ = L.cfield

    # Compute Wick contractions for each component
    keldysh = wick_contraction(ψ(Out()) * ψ'(In()), L, order; simplify, kwargs...)
    retarded = wick_contraction(ψ(Out()) * ϕ'(In()), L, order; simplify, kwargs...)
    advanced = wick_contraction(ϕ(Out()) * ψ'(In()), L, order; simplify, kwargs...)

    # Apply filtering and simplification to all components
    for component in (keldysh, retarded, advanced)
        filter_nonzero!(component)
        _simplify_prefactors!(component)
    end

    return DressedPropagator(keldysh, retarded, advanced, order, parameters(L)^order)
end

"get parameter of the interaction lagrangian"
parameters(d::DressedPropagator) = d.parameter

"""
$(DocStringExtensions.TYPEDEF)

A structure representing a collection of dressed propagator each involving a
different parameter of the expansion.

# Fields
$(DocStringExtensions.FIELDS)

"""
struct DressedPropagatorSum{GS}
    "The arguments of the dressed propagator sum. Each involving a different parameter."
    arguments::Dict{SymbolicUtils.BasicSymbolic{Number},GS}
    "The order of the dressed propagator in the perturbation series"
    order::Int64
end

SymbolicUtils.arguments(d::DressedPropagatorSum) = d.arguments
order(d::DressedPropagatorSum) = d.order
parameters(d::DressedPropagatorSum) = map(G -> G.parameter, arguments(d))

function DressedPropagator(
    Ls::LagrangianSum,
    order=1;
    simplify::Union{Bool,Vector{Bool}}=Bool[
        !should_regularise(L.lagrangian) for L in arguments(Ls)
    ],
    kwargs...,
)
    ϕ = first(arguments(Ls)).qfield
    ψ = first(arguments(Ls)).cfield

    simplify = isa(simplify, Bool) ? fill(simplify, length(Ls)) : simplify

    # Compute Wick contractions for each component
    keldysh_pairs = wick_contraction(ψ(Out()) * ψ'(In()), Ls, order; simplify, kwargs...)
    retarded_pairs = wick_contraction(ψ(Out()) * ϕ'(In()), Ls, order; simplify, kwargs...)
    advanced_pairs = wick_contraction(ϕ(Out()) * ψ'(In()), Ls, order; simplify, kwargs...)

    # Apply filtering and simplification to all components

    dict = Dict(
        begin
            for component in
                last.((keldysh_pairs[idx], retarded_pairs[idx], advanced_pairs[idx]))
                filter_nonzero!(component)
                _simplify_prefactors!(component)
            end
            first(keldysh_pairs[idx]) => DressedPropagator(
                last.((keldysh_pairs[idx], retarded_pairs[idx], advanced_pairs[idx]))...,
                order,
                first(keldysh_pairs[idx]),
            )
        end for idx in eachindex(keldysh_pairs)
    )

    return DressedPropagatorSum(dict, order)
end
