
###########################
#  Interaction Lagrangian
###########################

"""
$(TYPEDEF)

Represents an interaction Lagrangian

# Fields
$(FIELDS)

# Constructor
$(TYPEDSIGNATURES)

Constructs an InteractionLagrangian from a given [`QTerm`](@ref) expression.

# Requirements
The constructor enforces several constraints on the input expression and throws `AssertionError` if any of them are not met:
- Must be a bulk term ([`is_bulk`](@ref))
- Must be conserved ([`is_conserved`](@ref))
- Must be physical ([`is_physical`](@ref))
- Can only contain up to two different fields
- Fields must have opposite contours

"""
struct InteractionLagrangian{T}
    "The Lagrangian expression as a [`QTerm`](@ref)"
    lagrangian::T
    "The quantum field destruction operator"
    qfield::Destroy
    "The classical field destruction operator"
    cfield::Destroy
    "The position of the interaction Lagrangian"
    position::Position
    function InteractionLagrangian(expr::QTerm)
        fields = _extract_unique_fields(expr)
        contours = contour_integers(fields)

        _assert_lagrangian(expr, fields, contours)

        q_idx = findfirst(iszero, contours)
        c_idx = findfirst(isone, contours)
        return new{typeof(expr)}(
            expr, fields[q_idx], fields[c_idx], position(fields[q_idx])
        )
    end
end # Does not have to be type stable, as it is called only once
function _extract_unique_fields(expr)
    fields = allfields(expr)
    set_reg_to_zero!(fields)
    unique_fields = unique(fields)
    filter(is_annihilation, unique_fields)
end
function _assert_lagrangian(expr, fields, contours)
    @assert is_bulk(expr) "An interaction Lagrangian only accepts bulk terms"
    @assert is_conserved(expr) "An interaction Lagrangian only accepts conserved terms"
    @assert is_physical(expr) "An interaction Lagrangian only accepts physical terms"
    @assert length(fields) <= 2 "An interaction Lagrangian only accepts up to two different fields"
    @assert unique(contours) == contours "An interaction Lagrangian only accepts fields with opposite contours"
    return nothing
end  # TODO what if only quantum or only classical

"get position of the interaction lagrangian"
position(L::InteractionLagrangian) = L.position

"""
    is_conserved(a::QTerm)

Checks if an expression [`QTerm`])(@ref) is conserved. A conserved expression is one that has equal numbers of creation and annihilation operators.

See also: [`is_physical`](@ref)
"""
function is_conserved(args_nc_::Union{Tuple{<:QSym,<:QSym},Vector{QSym}})
    n_destroy = Bool[arg isa Destroy for arg in args_nc_]
    if length(args_nc_) == 0 || iszero(sum(n_destroy))
        return false
    else
        length(args_nc_) == 1
        return length(args_nc_) ÷ sum(n_destroy) == 2
    end
end
is_conserved(a::QMul) = is_conserved(a.args_nc)
is_conserved(a::QAdd) = all(is_conserved(arg) for arg in arguments(a))
is_conserved(a::QSym) = false

"""
    is_physical(a::QTerm)

Checks if an expression [`QTerm`])(@ref) is physical. A physical expression is one that if it has an `In` position field it also has an `Out` position field and vice versa ([`Position`])(@ref). Furthermore, `In` position field can only creation fields ([`Create`](@ref)) and `Out` position field can only have annihilation fields ([`Destroy`](@ref)).

See also: [`is_conserved`](@ref)
"""
function is_physical(args_nc_::Vector{<:QField})
    positions = [position(f) for f in args_nc_]
    # checks if a mul has both in-out in a lagrangian
    in_out = has_in(positions) ? (has_out(positions) ? true : false) : true
    out_in = has_out(positions) ? (has_in(positions) ? true : false) : true
    itr = map(is_physical, args_nc_)
    physical = all(itr) # individual field are physical
    return physical && in_out && out_in
end
@inline is_physical(a::QMul) = is_physical(a.args_nc)
is_physical(a::QAdd) = all(is_physical(arg) for arg in arguments(a))
is_physical(a::Destroy) = !is_in(a)
is_physical(a::Create) = !is_out(a)

function (L::InteractionLagrangian)(i::Int)
    new_lagrangian = set_position(L.lagrangian, Bulk(i))
    return InteractionLagrangian(new_lagrangian)
end
