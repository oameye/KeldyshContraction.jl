##########################################
#       dressed green's function
##########################################
"""
$(DocStringExtensions.TYPEDEF)

A structure representing dressed propagator in the Retarded-Advanced-Keldysh basis
([`PropagatorType`](@ref)).

# Fields
$(DocStringExtensions.FIELDS)
where it assumed that the fields are of type `Union{SymbolicUtils.Symbolic{<:Number}, Number}`.

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
    "Parameters of the perturbation series"
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
function DressedPropagator(L::InteractionLagrangian, order=1; kwargs...)
    ϕ = L.qfield
    ψ = L.cfield

    # Compute Wick contractions for each component
    keldysh = wick_contraction(ψ(Out()) * ψ'(In()), L, order; kwargs...)
    retarded = wick_contraction(ψ(Out()) * ϕ'(In()), L, order; kwargs...)
    advanced = wick_contraction(ϕ(Out()) * ψ'(In()), L, order; kwargs...)

    # Apply filtering and simplification to all components
    for component in (keldysh, retarded, advanced)
        filter_nonzero!(component)
        _simplify_prefactors!(component)
    end

    return DressedPropagator(keldysh, retarded, advanced, order, parameter(L))
end

"get parameter of the interaction lagrangian"
parameter(d::DressedPropagator) = d.parameter
