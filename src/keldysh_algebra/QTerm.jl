########################
#    Multiplication
########################

"""
$(DocStringExtensions.TYPEDEF)

Represent a multiplication involving quantum fields  of [`QSym`](@ref) types.

$(DocStringExtensions.FIELDS)
"""
struct QMul{T<:Number} <: QTerm
    "The commutative prefactor."
    arg_c::T
    "A vector containing all [`QSym`](@ref) types."
    args_nc::Vector{QSym}
    function QMul(arg_c::T, args_nc::Vector{<:QSym}) where {T<:Number}
        if isequal(arg_c, 0)
            return new{T}(0, QSym[])
        else
            return new{T}(arg_c, args_nc)
        end
    end
    QMul(args_nc::Vector{<:QSym}) = new{Int}(1, args_nc)
    QMul(s) = new{Int64}(1, [s])
    QMul{T}(s) where {T} = new{T}(T(1), [s])
    QMul() = new{Int64}(0, QSym[])
    QMul{T}() where {T} = new{T}(T(0), QSym[])
end

Base.promote_rule(::Type{QMul{S}}, ::Type{QMul{T}}) where {S,T} = QMul{promote_rule(S, T)}
function Base.convert(::Type{QMul{T}}, x::QMul{S}) where {T<:Number,S<:Number}
    return QMul(convert(T, x.arg_c), x.args_nc)
end
Base.length(a::QMul) = length(a.args_nc)

SymbolicUtils.operation(::QMul) = (*)
"""
    arguments(a::QMul)

Return the vector of the factors of [`QMul`](@ref).
"""
SymbolicUtils.arguments(a::QMul) = vcat(a.arg_c, a.args_nc)

function TermInterface.maketerm(::Type{<:QMul}, ::typeof(*), args, metadata)
    args_c = filter(x -> !(x isa QField), args)
    args_nc = filter(x -> x isa QField, args)
    arg_c = isempty(args_c) ? 1 : *(args_c...)
    # isempty(args_nc) && return arg_c
    return QMul(arg_c, args_nc)
end

TermInterface.metadata(a::QMul) = nothing

function Base.adjoint(q::QMul)
    args_nc = map(adjoint, q.args_nc)
    sort!(args_nc; by=position)
    sort!(args_nc; by=ladder)
    return QMul(conj(q.arg_c), args_nc)
end

Base.zero(q::QMul) = QMul(0.0, QSym[])
Base.iszero(q::QMul) = iszero(q.arg_c)

"""
    is_bulk(q::QTerm)

Checks if a term is in the bulk. A term is bulk if it has no `In` or `Out` position fields ([`Position`](@ref)).
"""
function is_bulk(q::QMul)
    args = q.args_nc
    bool = false
    for f in args
        if is_bulk(f)
            bool = true
        else
            return false
        end
    end
    return bool
end
allfields(q::QMul) = q.args_nc

########################
#       Addition
########################

"""
    QAdd <: QTerm

Represent an addition involving [`QField`](@ref) and other types.
"""
struct QAdd{T<:Number} <: QTerm
    arguments::Vector{QMul{T}}
    function QAdd(args::Vector{<:QMul})
        vs = promote(args...)
        new{typeof(first(vs)).parameters[1]}(collect(vs))
    end
    function QAdd(args::Vector{QMul{T}}) where {T<:Number}
        new{T}(args)
    end
    function QAdd(args::Vector{<:QSym})
        new{Int64}([QMul(s) for s in args])
    end
    QAdd{T}() where {T} = new{T}([QMul{T}()])
end
Base.length(a::QAdd) = length(a.arguments)
SymbolicUtils.operation(::QAdd) = (+)
"""
    arguments(a::QAdd)

Return the vector of the arguments of [`QAdd`](@ref).
"""
SymbolicUtils.arguments(a::QAdd) = a.arguments
TermInterface.maketerm(::Type{<:QAdd}, ::typeof(+), args, metadata) = QAdd(args)
TermInterface.metadata(a::QAdd) = nothing

Base.adjoint(q::QAdd) = QAdd(map(adjoint, arguments(q)))
function allfields(q::QAdd)
    vfields = map(allfields, arguments(q))
    reduce(vcat, vfields)
end

Base.promote_rule(::Type{QAdd{S}}, ::Type{QAdd{T}}) where {S,T} = QAdd{promote_rule(S, T)}
function Base.convert(::Type{QAdd{T}}, xadd::QAdd{S}) where {T<:Number,S<:Number}
    return QAdd([QMul(convert(T, x.arg_c), x.args_nc) for x in xadd.arguments])
end

#########################
#       isequal
#########################

function Base.isequal(a::QMul, b::QMul)
    isequal(a.arg_c, b.arg_c) || return false
    length(a.args_nc) == length(b.args_nc) || return false
    for (arg_a, arg_b) in zip(a.args_nc, b.args_nc)
        isequal(arg_a, arg_b) || return false
    end
    return true
end

function Base.isequal(a::QAdd, b::QAdd)
    length(arguments(a)) == length(arguments(b)) || return false
    for (arg_a, arg_b) in zip(arguments(a), arguments(b))
        isequal(arg_a, arg_b) || return false
    end
    return true
end

function Base.isequal(a::QAdd, b::QSym)
    args = arguments(a)
    if length(args) == 1
        return isequal(args[1], b)
    end
    return false
end
function Base.isequal(b::QSym, a::QAdd)
    args = arguments(a)
    if length(args) == 1
        return isequal(args[1], b)
    end
    return false
end

function Base.isequal(a::QMul, b::Int)
    args = a.args_nc
    if isempty(args)
        return isequal(a.arg_c, b)
    end
    if iszero(QMul) && iszero(b)
        return true
    end
    return false
end
function Base.isequal(b::Int, a::QMul)
    args = a.args_nc
    if isempty(args)
        return isequal(a.arg_c, b)
    end
    if iszero(QMul) && iszero(b)
        return true
    end
    return false
end
function Base.isequal(a::QMul, b::QSym)
    args = a.args_nc
    if length(args) == 1 && isone(a.arg_c)
        return isequal(args[1], b)
    end
    return false
end
function Base.isequal(b::QSym, a::QMul)
    args = a.args_nc
    if length(args) == 1 && isone(a.arg_c)
        return isequal(args[1], b)
    end
    return false
end

#########################
#       Position
#########################

is_bulk(q::QAdd) = all(is_bulk(q) for q in arguments(q))

function set_position_mul(a::QMul, p::Position)::QMul
    args = QSym[arg(p) for arg in a.args_nc]
    return QMul(a.arg_c, args)
end

function set_position(a::QAdd, p::Position)::QAdd
    args = QMul[set_position_mul(arg, p) for arg in arguments(a)]
    return QAdd(args)
end
