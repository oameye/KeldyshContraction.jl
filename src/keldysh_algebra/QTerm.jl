########################
#    Multiplication
########################

"""Concrete product of fields with coefficient type `C` and statistics `S`."""
struct QMul{C<:Number,S<:Statistics} <: QTerm
    arg_c::C
    args_nc::Vector{Field{S}}
end

@inline _apply_exchange_sign(c::C, sign::Int8) where {C<:Number} = sign == 1 ? c : -c

function QMul(arg_c::C, args_nc::Vector{Field{S}}) where {C<:Number,S<:Statistics}
    if iszero(arg_c)
        return QMul{C,S}(zero(C), Field{S}[])
    end
    args = copy(args_nc)
    sign = canonicalize_fields!(args)
    coeff = _apply_exchange_sign(arg_c, sign)
    return QMul{C,S}(coeff, args)
end

QMul(args_nc::Vector{Field{S}}) where {S<:Statistics} = QMul(1, args_nc)
QMul(s::Field{S}) where {S<:Statistics} = QMul(1, Field{S}[s])
QMul{C,S}() where {C<:Number,S<:Statistics} = QMul{C,S}(zero(C), Field{S}[])

Base.length(a::QMul) = length(a.args_nc)
Base.iszero(q::QMul) = iszero(q.arg_c)
Base.isone(q::QMul) = isone(q.arg_c) && isempty(q.args_nc)
Base.zero(q::QMul{C,S}) where {C,S} = QMul{C,S}(zero(C), Field{S}[])
Base.one(q::QMul{C,S}) where {C,S} = QMul{C,S}(one(C), Field{S}[])
Base.zero(::Type{QMul{C,S}}) where {C,S} = QMul{C,S}(zero(C), Field{S}[])
Base.one(::Type{QMul{C,S}}) where {C,S} = QMul{C,S}(one(C), Field{S}[])

function Base.promote_rule(
    ::Type{QMul{C1,S}}, ::Type{QMul{C2,S}}
) where {C1<:Number,C2<:Number,S<:Statistics}
    return QMul{promote_type(C1, C2),S}
end
function Base.convert(
    ::Type{QMul{C,S}}, q::QMul{D,S}
) where {C<:Number,D<:Number,S<:Statistics}
    return QMul{C,S}(convert(C, q.arg_c), copy(q.args_nc))
end

coefficient(q::QMul) = q.arg_c
fields(q::QMul) = q.args_nc
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
    ::Type{<:QMul}, ::typeof(*), args::Vector{Field{S}}, metadata
) where {S<:Statistics}
    return QMul(1, args)
end
function TermInterface.maketerm(
    ::Type{<:QMul},
    ::typeof(*),
    args::Vector{Union{C,Field{S}}},
    metadata,
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
    return QMul(coeff, fs)
end

function bar(q::QMul{C,S}) where {C,S}
    return QMul(conj(q.arg_c), Field{S}[bar(f) for f in q.args_nc])
end

function is_bulk(q::QMul)
    isempty(q.args_nc) && return true
    return all(is_bulk, q.args_nc)
end

########################
#       Addition
########################

"""Concrete sum of homogeneous `QMul{C,S}` terms."""
struct QAdd{C<:Number,S<:Statistics} <: QTerm
    arguments::Vector{QMul{C,S}}
end

QAdd(args::Vector{QMul{C,S}}) where {C<:Number,S<:Statistics} = QAdd{C,S}(args)
function QAdd(args::Vector{Field{S}}) where {S<:Statistics}
    return QAdd{Int,S}(QMul{Int,S}[QMul(1, Field{S}[f]) for f in args])
end
QAdd{C,S}() where {C<:Number,S<:Statistics} =
    QAdd{C,S}(QMul{C,S}[QMul{C,S}(zero(C), Field{S}[])])

Base.length(a::QAdd) = length(a.arguments)
Base.iszero(a::QAdd) = all(iszero, a.arguments)
Base.isone(a::QAdd) = length(a.arguments) == 1 && isone(first(a.arguments))
Base.zero(a::QAdd{C,S}) where {C,S} = QAdd{C,S}()
Base.one(a::QAdd{C,S}) where {C,S} =
    QAdd{C,S}(QMul{C,S}[QMul{C,S}(one(C), Field{S}[])])
Base.zero(::Type{QAdd{C,S}}) where {C,S} = QAdd{C,S}()
Base.one(::Type{QAdd{C,S}}) where {C,S} =
    QAdd{C,S}(QMul{C,S}[QMul{C,S}(one(C), Field{S}[])])

SymbolicUtils.operation(::QAdd) = (+)
SymbolicUtils.arguments(a::QAdd) = a.arguments
TermInterface.metadata(::QAdd) = nothing
function TermInterface.maketerm(
    ::Type{<:QAdd}, ::typeof(+), args::Vector{QMul{C,S}}, metadata
) where {C<:Number,S<:Statistics}
    return QAdd(args)
end

coefficient(q::QAdd) = map(coefficient, q.arguments)
terms(q::QAdd) = q.arguments
function allfields(q::QAdd{C,S}) where {C,S}
    out = Field{S}[]
    for term in q.arguments
        append!(out, term.args_nc)
    end
    return out
end
function bar(q::QAdd{C,S}) where {C,S}
    return QAdd(QMul{C,S}[bar(term) for term in q.arguments])
end

function Base.promote_rule(
    ::Type{QAdd{C1,S}}, ::Type{QAdd{C2,S}}
) where {C1<:Number,C2<:Number,S<:Statistics}
    return QAdd{promote_type(C1, C2),S}
end
function Base.convert(
    ::Type{QAdd{C,S}}, q::QAdd{D,S}
) where {C<:Number,D<:Number,S<:Statistics}
    return QAdd(QMul{C,S}[convert(QMul{C,S}, term) for term in q.arguments])
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

Base.:(==)(a::QMul, b::QMul) = isequal(a, b)
Base.:(==)(a::QAdd, b::QAdd) = isequal(a, b)

#########################
#       Position
#########################

is_bulk(q::QAdd) = all(is_bulk, q.arguments)

function set_position_mul(a::QMul{C,S}, p::Position) where {C<:Number,S<:Statistics}
    args = Field{S}[f(p) for f in a.args_nc]
    return QMul(a.arg_c, args)
end
function set_position(a::QAdd{C,S}, p::Position) where {C<:Number,S<:Statistics}
    return QAdd(QMul{C,S}[set_position_mul(arg, p) for arg in a.arguments])
end
