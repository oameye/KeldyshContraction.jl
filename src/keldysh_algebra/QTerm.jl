########################
#    Multiplication
########################

"""
$(DocStringExtensions.TYPEDEF)

A concrete product of path-integral fields: a coefficient of type `C` times an ordered
vector of [`Field`](@ref)s with statistics `S`. Storage is homogeneous, so a `QMul` is
recursively concrete and its element type never depends on runtime contents.

The canonical form is unique: fields are kept sorted, and a zero coefficient is stored with
an empty field vector so that algebraic zero has one representation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct QMul{C<:Number,S<:Statistics} <: QTerm
    "Scalar coefficient multiplying the product"
    arg_c::C
    "Canonically ordered noncommutative factors"
    args_nc::Vector{Field{S}}

    function QMul{C,S}(
        arg_c::C, args_nc::Vector{Field{S}}, ::Val{:sort}
    ) where {C<:Number,S<:Statistics}
        if iszero(arg_c)
            return new{C,S}(zero(C), Field{S}[])
        end
        sign = canonicalize_fields!(args_nc)
        coeff = _apply_exchange_sign(arg_c, sign)
        return new{C,S}(coeff, args_nc)
    end

    function QMul{C,S}(
        arg_c::C, args_nc::Vector{Field{S}}, ::Val{:presorted}
    ) where {C<:Number,S<:Statistics}
        if iszero(arg_c)
            return new{C,S}(zero(C), Field{S}[])
        end
        return new{C,S}(arg_c, args_nc)
    end
end

@inline _apply_exchange_sign(c::C, sign::Int8) where {C<:Number} = sign == 1 ? c : -c

"""
$(DocStringExtensions.SIGNATURES)

Build a `QMul` by sorting `args_nc` into canonical order in place. The caller must have
produced `args_nc` for this call alone: the vector is stored, not copied, and is mutated.
"""
@inline function _qmul_sorting(
    arg_c::C, args_nc::Vector{Field{S}}
) where {C<:Number,S<:Statistics}
    return QMul{C,S}(arg_c, args_nc, Val(:sort))
end

"""
$(DocStringExtensions.SIGNATURES)

Build a `QMul` from a vector that is already in canonical order, skipping the sort. The
vector is stored without copying, so it may be shared with another `QMul` only because
neither ever reorders it: every mutating path goes through [`_qmul_sorting`](@ref).
"""
@inline function _qmul_presorted(
    arg_c::C, args_nc::Vector{Field{S}}
) where {C<:Number,S<:Statistics}
    return QMul{C,S}(arg_c, args_nc, Val(:presorted))
end

# Public vector construction is defensive: callers keep ownership of `args_nc`.
function QMul{C,S}(arg_c::C, args_nc::Vector{Field{S}}) where {C<:Number,S<:Statistics}
    return _qmul_sorting(arg_c, copy(args_nc))
end
function QMul(arg_c::C, args_nc::Vector{Field{S}}) where {C<:Number,S<:Statistics}
    return QMul{C,S}(arg_c, args_nc)
end

QMul(args_nc::Vector{Field{S}}) where {S<:Statistics} = QMul(1, args_nc)
QMul(s::Field{S}) where {S<:Statistics} = _qmul_sorting(1, Field{S}[s])
QMul{C,S}() where {C<:Number,S<:Statistics} = _qmul_sorting(zero(C), Field{S}[])

Base.length(a::QMul) = length(a.args_nc)
Base.iszero(q::QMul) = iszero(q.arg_c)
Base.isone(q::QMul) = isone(q.arg_c) && isempty(q.args_nc)
Base.zero(q::QMul{C,S}) where {C,S} = _qmul_sorting(zero(C), Field{S}[])
Base.one(q::QMul{C,S}) where {C,S} = _qmul_sorting(one(C), Field{S}[])
Base.zero(::Type{QMul{C,S}}) where {C,S} = _qmul_sorting(zero(C), Field{S}[])
Base.one(::Type{QMul{C,S}}) where {C,S} = _qmul_sorting(one(C), Field{S}[])

function Base.promote_rule(
    ::Type{QMul{C1,S}}, ::Type{QMul{C2,S}}
) where {C1<:Number,C2<:Number,S<:Statistics}
    return QMul{promote_type(C1, C2),S}
end
Base.convert(::Type{QMul{C,S}}, q::QMul{C,S}) where {C<:Number,S<:Statistics} = q
function Base.convert(
    ::Type{QMul{C,S}}, q::QMul{D,S}
) where {C<:Number,D<:Number,S<:Statistics}
    return _qmul_presorted(convert(C, q.arg_c), q.args_nc)
end

"""Scalar coefficient of a product."""
coefficient(q::QMul) = q.arg_c

"""
Canonically ordered fields of a product, as a copy the caller may mutate freely. Internal
call sites that must not allocate read `q.args_nc` directly.
"""
fields(q::QMul) = copy(q.args_nc)

"""Products making up an expression; a `QMul` is a one-term sum."""
terms(q::QMul) = (q,)

allfields(q::QMul) = copy(q.args_nc)

SymbolicUtils.operation(::QMul) = (*)
function SymbolicUtils.arguments(a::QMul{C,S}) where {C,S}
    out = Vector{Union{C,Field{S}}}(undef, length(a.args_nc) + 1)
    out[1] = a.arg_c
    for i in eachindex(a.args_nc)
        out[i + 1] = a.args_nc[i]
    end
    return out
end
TermInterface.metadata(::QMul) = nothing

function TermInterface.maketerm(
    ::Type{QMul{C,S}}, ::typeof(*), args::Vector{Field{S}}, metadata
) where {C<:Number,S<:Statistics}
    return QMul{C,S}(one(C), args)
end
function TermInterface.maketerm(
    ::Type{QMul{C,S}}, ::typeof(*), args::Vector{Union{C,Field{S}}}, metadata
) where {C<:Number,S<:Statistics}
    coeff = one(C)
    fs = Field{S}[]
    for arg in args
        if arg isa Field{S}
            push!(fs, arg)
        else
            coeff *= arg
        end
    end
    return _qmul_sorting(coeff, fs)
end

function bar(q::QMul{C,S}) where {C,S}
    return _qmul_sorting(conj(q.arg_c), Field{S}[bar(f) for f in q.args_nc])
end

function is_bulk(q::QMul)
    isempty(q.args_nc) && return true
    return all(is_bulk, q.args_nc)
end

########################
#       Addition
########################

"""
$(DocStringExtensions.TYPEDEF)

A concrete sum of homogeneous [`QMul`](@ref) terms sharing one coefficient type `C` and one
statistics `S`. Storage is a `Vector{QMul{C,S}}`, so a `QAdd` is recursively concrete.

A `QAdd` always holds at least one term: algebraic zero is the single zero-coefficient
monomial, never an empty vector. Hence `iszero`, `isone`, and `first(terms(q))` are total.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct QAdd{C<:Number,S<:Statistics} <: QTerm
    "Summands, at least one"
    arguments::Vector{QMul{C,S}}

    function QAdd{C,S}(
        arguments::Vector{QMul{C,S}}, ::Val{:owned}
    ) where {C<:Number,S<:Statistics}
        isempty(arguments) && throw(
            ArgumentError(
                "a QAdd needs at least one term; algebraic zero is zero(QAdd{$C,$S})"
            ),
        )
        return new{C,S}(arguments)
    end
end

"""Construct a `QAdd` from an internally owned term vector without copying it."""
@inline function _qadd_owned(arguments::Vector{QMul{C,S}}) where {C<:Number,S<:Statistics}
    return QAdd{C,S}(arguments, Val(:owned))
end

# Public vector construction is defensive: callers retain ownership of `arguments`.
function QAdd{C,S}(arguments::Vector{QMul{C,S}}) where {C<:Number,S<:Statistics}
    return _qadd_owned(copy(arguments))
end
QAdd(arguments::Vector{QMul{C,S}}) where {C<:Number,S<:Statistics} = QAdd{C,S}(arguments)

function QAdd(args::Vector{Field{S}}) where {S<:Statistics}
    return _qadd_owned(QMul{Int,S}[_qmul_sorting(1, Field{S}[f]) for f in args])
end
function QAdd{C,S}() where {C<:Number,S<:Statistics}
    return _qadd_owned(QMul{C,S}[_qmul_sorting(zero(C), Field{S}[])])
end

Base.length(a::QAdd) = length(a.arguments)
Base.iszero(a::QAdd) = all(iszero, a.arguments)
Base.isone(a::QAdd) = length(a.arguments) == 1 && isone(first(a.arguments))
Base.zero(a::QAdd{C,S}) where {C,S} = QAdd{C,S}()
Base.one(a::QAdd{C,S}) where {C,S} =
    _qadd_owned(QMul{C,S}[_qmul_sorting(one(C), Field{S}[])])
Base.zero(::Type{QAdd{C,S}}) where {C,S} = QAdd{C,S}()
function Base.one(::Type{QAdd{C,S}}) where {C,S}
    return _qadd_owned(QMul{C,S}[_qmul_sorting(one(C), Field{S}[])])
end

SymbolicUtils.operation(::QAdd) = (+)
SymbolicUtils.arguments(a::QAdd) = copy(a.arguments)
TermInterface.metadata(::QAdd) = nothing
function TermInterface.maketerm(
    ::Type{QAdd{C,S}}, ::typeof(+), args::Vector{QMul{C,S}}, metadata
) where {C<:Number,S<:Statistics}
    return QAdd{C,S}(args)
end

"""Coefficients of every term of a sum."""
coefficient(q::QAdd) = map(coefficient, q.arguments)

"""Products making up a sum, as a copy the caller may mutate freely."""
terms(q::QAdd) = copy(q.arguments)

"""Fields appearing anywhere in a sum, in term order, as a fresh vector."""
fields(q::QAdd) = allfields(q)
function allfields(q::QAdd{C,S}) where {C,S}
    out = Field{S}[]
    for term in q.arguments
        append!(out, term.args_nc)
    end
    return out
end
function bar(q::QAdd{C,S}) where {C,S}
    return _qadd_owned(QMul{C,S}[bar(term) for term in q.arguments])
end

function Base.promote_rule(
    ::Type{QAdd{C1,S}}, ::Type{QAdd{C2,S}}
) where {C1<:Number,C2<:Number,S<:Statistics}
    return QAdd{promote_type(C1, C2),S}
end
Base.convert(::Type{QAdd{C,S}}, q::QAdd{C,S}) where {C<:Number,S<:Statistics} = q
function Base.convert(
    ::Type{QAdd{C,S}}, q::QAdd{D,S}
) where {C<:Number,D<:Number,S<:Statistics}
    return _qadd_owned(QMul{C,S}[convert(QMul{C,S}, term) for term in q.arguments])
end

#########################
#       isequal
#########################

function Base.isequal(a::QMul, b::QMul)
    isequal(a.arg_c, b.arg_c) || return false
    length(a.args_nc) == length(b.args_nc) || return false
    return all(isequal(x, y) for (x, y) in zip(a.args_nc, b.args_nc))
end

function Base.isequal(a::QAdd, b::QAdd)
    length(a.arguments) == length(b.arguments) || return false
    return all(isequal(x, y) for (x, y) in zip(a.arguments, b.arguments))
end

function Base.isequal(a::QAdd, b::Field)
    length(a.arguments) == 1 || return false
    return isequal(first(a.arguments), b)
end
Base.isequal(b::Field, a::QAdd) = isequal(a, b)

function Base.isequal(a::QMul, b::Number)
    return isempty(a.args_nc) && isequal(a.arg_c, b)
end
Base.isequal(b::Number, a::QMul) = isequal(a, b)

function Base.isequal(a::QMul, b::Field)
    return length(a.args_nc) == 1 && isone(a.arg_c) && isequal(first(a.args_nc), b)
end
Base.isequal(b::Field, a::QMul) = isequal(a, b)

# `==` mirrors `isequal` for every pair `isequal` accepts. Defining only the homogeneous
# pairs made `one(f) == 1` silently false while `isequal(one(f), 1)` was true.
Base.:(==)(a::QMul, b::QMul) = isequal(a, b)
Base.:(==)(a::QAdd, b::QAdd) = isequal(a, b)
Base.:(==)(a::QMul, b::Number) = isequal(a, b)
Base.:(==)(a::Number, b::QMul) = isequal(a, b)
Base.:(==)(a::QMul, b::Field) = isequal(a, b)
Base.:(==)(a::Field, b::QMul) = isequal(a, b)
Base.:(==)(a::QAdd, b::Field) = isequal(a, b)
Base.:(==)(a::Field, b::QAdd) = isequal(a, b)

#########################
#       Position
#########################

is_bulk(q::QAdd) = all(is_bulk, q.arguments)

function set_position_mul(a::QMul{C,S}, p::Position) where {C<:Number,S<:Statistics}
    args = Field{S}[f(p) for f in a.args_nc]
    return _qmul_sorting(a.arg_c, args)
end
function set_position(a::QAdd{C,S}, p::Position) where {C<:Number,S<:Statistics}
    return _qadd_owned(QMul{C,S}[set_position_mul(arg, p) for arg in a.arguments])
end
