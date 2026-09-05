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

Represents an interaction Lagrangian.

# Fields
$(DocStringExtensions.FIELDS)

# Constructor
$(DocStringExtensions.TYPEDSIGNATURES)

Constructs an `InteractionLagrangian` from a [`QTerm`](@ref) expression.

# Requirements
The expression must:
- be a bulk term (`is_bulk`);
- be conserved ([`is_conserved`](@ref));
- be physical ([`is_physical`](@ref));
- contain both Keldysh components.
"""
struct InteractionLagrangian{C<:Number,S<:Statistics} <: Lagrangian
    "The Lagrangian expression"
    lagrangian::QAdd{C,S}
    "Physical field families appearing in the interaction"
    families::Vector{FieldFamily{S}}
    "The position of the interaction Lagrangian"
    position::Position
    "Parameter of the perturbation series"
    parameter::CSym
end

normalize_interaction(expr::QAdd) = expr
function normalize_interaction(expr::QMul{C,S}) where {C<:Number,S<:Statistics}
    return QAdd(QMul{C,S}[expr])
end

function InteractionLagrangian(
    expr::Union{QMul{C,S},QAdd{C,S}}, parameter=DEFAULT_PARAMETER
) where {C<:Number,S<:Statistics}
    expression = normalize_interaction(expr)
    fields = allfields(expression)
    assert_lagrangian(expression, fields)

    return InteractionLagrangian{C,S}(
        expression,
        interaction_families(expression),
        position(first(fields)),
        parameter,
    )
end

function interaction_families(expr::Union{QMul{C,S},QAdd{C,S}}) where {C<:Number,S<:Statistics}
    families = unique(field_family.(allfields(expr)))
    sort!(families)
    return families
end

function assert_lagrangian(expr, fields)
    is_bulk(expr) ||
        throw(ArgumentError("an interaction Lagrangian only accepts bulk terms"))
    is_conserved(expr) ||
        throw(ArgumentError("an interaction Lagrangian only accepts conserved terms"))
    is_physical(expr) ||
        throw(ArgumentError("an interaction Lagrangian only accepts physical terms"))
    contours = contour_integers(fields)
    any(iszero, contours) && any(isone, contours) || throw(
        ArgumentError("an interaction Lagrangian must contain both Keldysh components")
    )
    return nothing
end

position(L::InteractionLagrangian) = L.position
parameters(L::InteractionLagrangian) = L.parameter
field_families(L::InteractionLagrangian) = copy(L.families)

function target_family(L::InteractionLagrangian)
    length(L.families) == 1 || throw(
        ArgumentError("an interaction with multiple field families requires an explicit target")
    )
    return only(L.families)
end

function target_family(L::InteractionLagrangian{C,S}, target::FieldFamily{S}) where {C,S}
    target in L.families || throw(ArgumentError("target field family is not in the interaction"))
    return target
end

# Keep tuple/vector methods separate so the statistics parameter remains bound.
function balanced_orientation(args)
    n_unbarred = count(is_unbarred, args)
    return n_unbarred > 0 && n_unbarred == count(is_barred, args)
end

"""
    is_conserved(a)

Checks whether a field expression is conserved. A conserved expression contains equal,
nonzero numbers of barred and unbarred fields.

See also: [`is_physical`](@ref)
"""
is_conserved(args::Vector{Field{S}}) where {S<:Statistics} = balanced_orientation(args)
is_conserved(::Tuple{}) = false
function is_conserved(args::Tuple{Field{S},Vararg{Field{S}}}) where {S<:Statistics}
    return balanced_orientation(args)
end
is_conserved(a::QMul) = is_conserved(a.args_nc)
is_conserved(a::QAdd) = all(is_conserved, a.arguments)
is_conserved(::Field) = false

"""
    is_physical(a)

Checks whether a field expression is physical. If an expression contains an [`In`](@ref)
field, it must also contain an [`Out`](@ref) field and vice versa. An `In` field must be
barred, while an `Out` field must be unbarred.

See also: [`is_conserved`](@ref)
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

Represents a sum of interaction Lagrangians, each describing a different process.
All terms must use the same physical field families.

# Fields
$(DocStringExtensions.FIELDS)

# Constructor
$(DocStringExtensions.TYPEDSIGNATURES)
"""
struct LagrangianSum{C<:Number,S<:Statistics} <: Lagrangian
    "Interaction Lagrangians in the sum"
    arguments::Vector{InteractionLagrangian{C,S}}

    function LagrangianSum(
        arguments::Vector{InteractionLagrangian{C,S}}
    ) where {C<:Number,S<:Statistics}
        isempty(arguments) && throw(ArgumentError("a LagrangianSum cannot be empty"))
        check_common_fields(arguments)
        return new{C,S}(arguments)
    end
end

function check_common_fields(args::AbstractVector{<:InteractionLagrangian})
    allequal(getfield.(args, :families)) || throw(
        ArgumentError("all InteractionLagrangian must have the same physical field families")
    )
    return nothing
end

field_families(Ls::LagrangianSum) = field_families(first(arguments(Ls)))
target_family(Ls::LagrangianSum) = target_family(first(arguments(Ls)))
function target_family(Ls::LagrangianSum{C,S}, target::FieldFamily{S}) where {C,S}
    return target_family(first(arguments(Ls)), target)
end

SymbolicUtils.arguments(Ls::LagrangianSum) = Ls.arguments
parameters(Ls::LagrangianSum) = [L.parameter for L in arguments(Ls)]
Base.length(a::LagrangianSum) = length(arguments(a))
SymbolicUtils.operation(::LagrangianSum) = (+)

function Base.:+(
    a::InteractionLagrangian{C1,S}, b::InteractionLagrangian{C2,S}
) where {C1<:Number,C2<:Number,S<:Statistics}
    C = promote_type(C1, C2)
    return LagrangianSum(
        InteractionLagrangian{C,S}[
            convert(InteractionLagrangian{C,S}, a),
            convert(InteractionLagrangian{C,S}, b),
        ]
    )
end

function Base.promote_rule(
    ::Type{InteractionLagrangian{C1,S}}, ::Type{InteractionLagrangian{C2,S}}
) where {C1<:Number,C2<:Number,S<:Statistics}
    return InteractionLagrangian{promote_type(C1, C2),S}
end

function Base.convert(
    ::Type{InteractionLagrangian{C,S}}, L::InteractionLagrangian{C,S}
) where {C<:Number,S<:Statistics}
    return L
end
function Base.convert(
    ::Type{InteractionLagrangian{C,S}}, L::InteractionLagrangian{D,S}
) where {C<:Number,D<:Number,S<:Statistics}
    return InteractionLagrangian{C,S}(
        convert(QAdd{C,S}, L.lagrangian), copy(L.families), L.position, L.parameter
    )
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
