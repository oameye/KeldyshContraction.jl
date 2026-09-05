const DEFAULT_PARAMETER = let
    @syms g::Number
    g
end

abstract type Lagrangian end

###########################
#  Interaction Lagrangian
###########################

"""Represents a bosonic interaction Lagrangian."""
struct InteractionLagrangian{T} <: Lagrangian
    lagrangian::T
    qfield::Field{Boson}
    cfield::Field{Boson}
    position::Position
    parameter::CSym
end

function InteractionLagrangian(
    expr::Union{QMul{C,Boson},QAdd{C,Boson}}, parameter=DEFAULT_PARAMETER
) where {C<:Number}
    fields = _extract_unique_fields(expr)
    contours = contour_integers(fields)

    _assert_lagrangian(expr, fields, contours)

    q_idx = findfirst(iszero, contours)
    c_idx = findfirst(isone, contours)
    q_idx === nothing && throw(ArgumentError("interaction has no quantum field"))
    c_idx === nothing && throw(ArgumentError("interaction has no classical field"))

    return InteractionLagrangian{typeof(expr)}(
        expr, fields[q_idx], fields[c_idx], position(fields[q_idx]), parameter
    )
end

function _extract_unique_fields(expr::Union{QMul{C,Boson},QAdd{C,Boson}}) where {C<:Number}
    fs = allfields(expr)
    set_reg_to_zero!(fs)
    unique_fields = unique(fs)
    result = Field{Boson}[]
    sizehint!(result, length(unique_fields))
    for field in unique_fields
        is_unbarred(field) && push!(result, field)
    end
    return result
end

function _assert_lagrangian(expr, fields, contours)
    is_bulk(expr) ||
        throw(ArgumentError("an interaction Lagrangian only accepts bulk terms"))
    is_conserved(expr) ||
        throw(ArgumentError("an interaction Lagrangian only accepts conserved terms"))
    is_physical(expr) ||
        throw(ArgumentError("an interaction Lagrangian only accepts physical terms"))
    length(fields) <= 2 || throw(
        ArgumentError("an interaction Lagrangian only accepts up to two different fields"),
    )
    unique(contours) == contours || throw(
        ArgumentError(
            "an interaction Lagrangian only accepts fields with opposite contours"
        ),
    )
    return nothing
end

position(L::InteractionLagrangian) = L.position
parameters(L::InteractionLagrangian) = L.parameter

# Keep tuple/vector methods separate so the statistics parameter remains bound.
function balanced_orientation(args)
    n_unbarred = count(is_unbarred, args)
    return n_unbarred > 0 && n_unbarred == count(is_barred, args)
end

"""Check whether fields contain equal, nonzero barred and unbarred counts."""
is_conserved(args::Vector{Field{S}}) where {S<:Statistics} = balanced_orientation(args)
is_conserved(::Tuple{}) = false
function is_conserved(args::Tuple{Field{S},Vararg{Field{S}}}) where {S<:Statistics}
    return balanced_orientation(args)
end
is_conserved(a::QMul) = is_conserved(a.args_nc)
is_conserved(a::QAdd) = all(is_conserved, a.arguments)
is_conserved(::Field) = false

"""Check whether fields satisfy the external endpoint convention."""
function is_physical(args::Vector{Field{S}}) where {S<:Statistics}
    positions = Position[position(f) for f in args]
    in_out = !has_in(positions) || has_out(positions)
    out_in = !has_out(positions) || has_in(positions)
    return all(is_physical, args) && in_out && out_in
end
is_physical(a::QMul) = is_physical(a.args_nc)
is_physical(a::QAdd) = all(is_physical, a.arguments)
is_physical(f::Field) = is_unbarred(f) ? !is_in(f) : !is_out(f)

function (L::InteractionLagrangian)(i::Int)
    new_lagrangian = set_position(L.lagrangian, Bulk(i))
    return InteractionLagrangian(new_lagrangian, parameters(L))
end

###########################
#      LagrangianSum
###########################

"""Sum of interaction Lagrangians with common fields."""
struct LagrangianSum{T} <: Lagrangian
    arguments::Vector{InteractionLagrangian{T}}
    function LagrangianSum(args::Vector{InteractionLagrangian{T}}, ::Val{:raw}) where {T}
        return new{T}(args)
    end
end

@unstable function LagrangianSum(args::Vector{<:InteractionLagrangian})
    check_common_fields(args)
    vs = promote(args...)
    return LagrangianSum(collect(vs))
end
function LagrangianSum(args::Vector{InteractionLagrangian{T}}) where {T}
    check_common_fields(args)
    return LagrangianSum(args, Val(:raw))
end

function check_common_fields(args::AbstractVector{<:InteractionLagrangian})
    allequal(getfield.(args, :qfield)) ||
        throw(ArgumentError("all InteractionLagrangian must have the same quantum field"))
    allequal(getfield.(args, :cfield)) ||
        throw(ArgumentError("all InteractionLagrangian must have the same classical field"))
    return nothing
end

SymbolicUtils.arguments(Ls::LagrangianSum) = Ls.arguments
parameters(Ls::LagrangianSum) = [L.parameter for L in arguments(Ls)]
Base.length(a::LagrangianSum) = length(arguments(a))
SymbolicUtils.operation(::LagrangianSum) = (+)

function Base.:+(a::InteractionLagrangian{T}, b::InteractionLagrangian{S}) where {T,S}
    P = promote_type(T, S)
    args = InteractionLagrangian{P}[convert(InteractionLagrangian{P}, a)]
    push!(args, convert(InteractionLagrangian{P}, b))
    return LagrangianSum(args)
end

function Base.promote_rule(
    ::Type{InteractionLagrangian{S}}, ::Type{InteractionLagrangian{T}}
) where {S,T}
    return InteractionLagrangian{promote_type(S, T)}
end
function Base.convert(
    ::Type{InteractionLagrangian{T}}, L::InteractionLagrangian{S}
) where {T,S}
    return InteractionLagrangian(convert(T, L.lagrangian), parameters(L))
end

#############################
# Explicit coefficient APIs
#############################

"""Convert all numeric coefficients in `q` to type `C`."""
function convert_coefficients(
    ::Type{C}, q::QMul{D,S}
) where {C<:Number,D<:Number,S<:Statistics}
    return QMul(convert(C, q.arg_c), copy(q.args_nc))
end
function convert_coefficients(
    ::Type{C}, q::QAdd{D,S}
) where {C<:Number,D<:Number,S<:Statistics}
    return QAdd(QMul{C,S}[convert_coefficients(C, term) for term in q.arguments])
end

"""Rationalize floating-point coefficients in `q`."""
rationalize_coefficients(q::QMul{C,S}) where {C<:Integer,S<:Statistics} = q
rationalize_coefficients(q::QMul{C,S}) where {C<:Rational,S<:Statistics} = q
rationalize_coefficients(q::QMul{Complex{C},S}) where {C<:Rational,S<:Statistics} = q

function rationalize_coefficients(q::QMul{C,S}) where {C<:AbstractFloat,S<:Statistics}
    R = typeof(rationalize(zero(C)))
    return QMul{R,S}(rationalize(q.arg_c), copy(q.args_nc))
end
function rationalize_coefficients(
    q::QMul{Complex{C},S}
) where {C<:AbstractFloat,S<:Statistics}
    R = typeof(rationalize(zero(C)))
    CR = Complex{R}
    c = complex(rationalize(real(q.arg_c)), rationalize(imag(q.arg_c)))
    return QMul{CR,S}(convert(CR, c), copy(q.args_nc))
end
rationalize_coefficients(q::QMul{C,S}) where {C<:Number,S<:Statistics} = q

function rationalize_coefficients(q::QAdd{C,S}) where {C<:Number,S<:Statistics}
    converted = map(rationalize_coefficients, q.arguments)
    return QAdd(converted)
end
