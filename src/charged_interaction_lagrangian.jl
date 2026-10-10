"""
$(DocStringExtensions.TYPEDEF)

Interaction Lagrangian with an explicit integer U(1) charge assignment for each field family.

`ChargedInteractionLagrangian` is intended for interactions such as a Hubbard--Stratonovich
cubic vertex, where every monomial is neutral under a weighted charge although it need not
contain equal raw numbers of barred and unbarred fields. The ordinary
[`InteractionLagrangian`](@ref) contract is unchanged.

Construct with charge pairs and an optional perturbative parameter, for example

```julia
ChargedInteractionLagrangian(expr, ψ => 1, χ => 2; parameter=:h)
```

A barred field contributes the negative of its family's charge and an unbarred field the positive
charge. Every monomial must have total charge zero.
"""
struct ChargedInteractionLagrangian{C<:Number,S<:Statistics,N} <: Lagrangian
    lagrangian::QAdd{C,S}
    families::Vector{FieldFamily{S}}
    charges::NTuple{N,Pair{FieldFamily{S},Int}}
    position::Position
    parameter::ParameterMonomial
end

function _canonical_charge_assignment(
    families::Vector{FieldFamily{S}}, assignments::Tuple, ::Val{N}
) where {S<:Statistics,N}
    N > 0 || throw(ArgumentError("a charged interaction requires field charges"))

    charges = Pair{FieldFamily{S},Int}[]
    seen = Set{FieldFamily{S}}()
    for assignment in assignments
        assignment isa Pair || throw(
            ArgumentError("field charges must be supplied as family => integer pairs")
        )
        family, charge = assignment
        family isa FieldFamily{S} || throw(
            ArgumentError("charge assignment statistics do not match the interaction")
        )
        charge isa Integer || throw(ArgumentError("field charges must be integers"))
        family in seen && throw(
            ArgumentError("duplicate charge assignment for field family $(name(family))"),
        )
        push!(seen, family)
        push!(charges, family => Int(charge))
    end

    length(charges) == length(families) || throw(
        ArgumentError(
            "every interaction field family requires exactly one charge assignment"
        ),
    )
    all(family -> family in seen, families) || throw(
        ArgumentError(
            "charge assignments must match the interaction field families exactly"
        ),
    )
    sort!(charges; by=first)
    return ntuple(index -> charges[index], Val(N))
end

function _assigned_charge(
    charges::NTuple{N,Pair{FieldFamily{S},Int}}, family::FieldFamily{S}
) where {S,N}
    for assignment in charges
        isequal(first(assignment), family) && return last(assignment)
    end
    return throw(ArgumentError("field family has no charge assignment"))
end

function _weighted_charge(args, charges)
    total = 0
    for field in args
        charge = _assigned_charge(charges, field_family(field))
        total += is_unbarred(field) ? charge : -charge
    end
    return total
end

function _assert_charged_lagrangian(expr, fields, charges)
    is_bulk(expr) ||
        throw(ArgumentError("a charged interaction Lagrangian only accepts bulk terms"))
    is_physical(expr) ||
        throw(ArgumentError("a charged interaction Lagrangian only accepts physical terms"))
    all(term -> iszero(_weighted_charge(term.args_nc, charges)), terms(expr)) || throw(
        ArgumentError(
            "every charged interaction monomial must have zero weighted U(1) charge"
        ),
    )
    contours = contour_integers(fields)
    any(iszero, contours) && any(isone, contours) || throw(
        ArgumentError(
            "a charged interaction Lagrangian must contain both Keldysh components"
        ),
    )
    return nothing
end

function ChargedInteractionLagrangian(
    expr::Union{QMul{C,S},QAdd{C,S}},
    assignments::Vararg{Pair,N};
    parameter=DEFAULT_PARAMETER,
) where {C<:Number,S<:Statistics,N}
    expression = normalize_interaction(expr)
    fields = allfields(expression)
    families = interaction_families(expression)
    charges = _canonical_charge_assignment(families, assignments, Val(N))
    _assert_charged_lagrangian(expression, fields, charges)

    return ChargedInteractionLagrangian{C,S,N}(
        expression,
        families,
        charges,
        position(first(fields)),
        parameter_monomial(parameter),
    )
end

position(L::ChargedInteractionLagrangian) = L.position
parameters(L::ChargedInteractionLagrangian) = L.parameter
field_families(L::ChargedInteractionLagrangian) = copy(L.families)

"""Return the integer U(1) charge assigned to `family` in a charged interaction."""
function field_charge(
    L::ChargedInteractionLagrangian{C,S,N}, family::FieldFamily{S}
) where {C,S,N}
    family in L.families ||
        throw(ArgumentError("field family is not in the charged interaction"))
    return _assigned_charge(L.charges, family)
end

function target_family(L::ChargedInteractionLagrangian)
    length(L.families) == 1 || throw(
        ArgumentError(
            "a charged interaction with multiple field families requires an explicit target",
        ),
    )
    return only(L.families)
end

function target_family(
    L::ChargedInteractionLagrangian{C,S,N}, target::FieldFamily{S}
) where {C,S,N}
    target in L.families ||
        throw(ArgumentError("target field family is not in the charged interaction"))
    return target
end

function (L::ChargedInteractionLagrangian)(i::Int)
    new_lagrangian = set_position(L.lagrangian, Bulk(i))
    return ChargedInteractionLagrangian(
        new_lagrangian, L.charges...; parameter=parameters(L)
    )
end
