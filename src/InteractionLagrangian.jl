const DEFAULT_PARAMETER = let
    @syms g::Number
    g
end

abstract type Lagrangian end

###########################
#  Interaction Lagrangian
###########################

"""
$(DocStringExtensions.TYPEDEF)

Represents a bosonic interaction Lagrangian.

The representation is migrated to `Field{Boson}` here; the field-family-neutral storage
redesign belongs to #241. Construction preserves the coefficient representation of the
input expression: exact/rational conversion is always explicit, via
[`rationalize_coefficients`](@ref) or [`convert_coefficients`](@ref).

# Fields
$(DocStringExtensions.FIELDS)
"""
struct InteractionLagrangian{T} <: Lagrangian
    "The symbolic interaction expression"
    lagrangian::T
    "The quantum field of the interaction"
    qfield::Field{Boson}
    "The classical field of the interaction"
    cfield::Field{Boson}
    "Bulk coordinate at which the interaction sits"
    position::Position
    "Symbolic coupling parameter"
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

# Shared body. The public methods below stay separate rather than uniting behind one
# `Union`, so that `S` stays bound: an empty `Tuple{Vararg{Field{S}}}` satisfies any `S`
# and would trip `Aqua.test_unbound_args`.
function _is_conserved(args)
    n_unbarred = count(is_unbarred, args)
    return n_unbarred > 0 && n_unbarred == count(is_barred, args)
end

"""
Return whether a collection of fields contains equal, nonzero counts of barred and unbarred
fields. An empty collection is not conserved: there is nothing to conserve, and a single
field is never conserved on its own.
"""
is_conserved(args::Vector{Field{S}}) where {S<:Statistics} = _is_conserved(args)
is_conserved(::Tuple{}) = false
function is_conserved(args::Tuple{Field{S},Vararg{Field{S}}}) where {S<:Statistics}
    return _is_conserved(args)
end
is_conserved(a::QMul) = is_conserved(a.args_nc)
is_conserved(a::QAdd) = all(is_conserved, a.arguments)
is_conserved(::Field) = false

"""
Check the current physical endpoint convention: an `In` field is barred and an `Out`
field is unbarred. Bulk fields of either orientation are physical.
"""
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

"""
$(DocStringExtensions.TYPEDEF)

A sum of [`InteractionLagrangian`](@ref)s sharing one quantum and one classical field.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LagrangianSum{T} <: Lagrangian
    "The summed interaction Lagrangians"
    arguments::Vector{InteractionLagrangian{T}}
    function LagrangianSum(args::Vector{InteractionLagrangian{T}}, ::Val{:raw}) where {T}
        return new{T}(args)
    end
end

# This constructor is migrated fully in #241. It remains isolated from the symbolic IR
# and is one of the pre-existing package-wide instability sites tracked by #238.
@unstable function LagrangianSum(args::Vector{<:InteractionLagrangian})
    _assert_common_fields(args)
    vs = promote(args...)
    return LagrangianSum(collect(vs))
end
function LagrangianSum(args::Vector{InteractionLagrangian{T}}) where {T}
    _assert_common_fields(args)
    return LagrangianSum(args, Val(:raw))
end

"""Every summand must share one quantum and one classical field."""
function _assert_common_fields(args::AbstractVector{<:InteractionLagrangian})
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

"""
    convert_coefficients(C, q)

Convert every numeric coefficient in a symbolic field expression `q` to the numeric type
`C`, while preserving its fields and statistics.
"""
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

"""
    rationalize_coefficients(q)

Explicitly replace floating-point coefficients in `q` by their `rationalize`d exact
representations. Exact integer and rational coefficients are returned unchanged.
"""
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
