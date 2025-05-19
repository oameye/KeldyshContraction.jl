########################
#    Multiplication
########################

"""
$(TYPEDEF)

Represent a multiplication involving quantum fields  of [`QSym`](@ref) types.

$(FIELDS)
"""
struct QMul{T<:Number} <: QTerm
    "The commutative prefactor."
    arg_c::T
    "A vector containing all [`QSym`](@ref) types."
    args_nc::Vector{QSym}
    function QMul(arg_c::T, args_nc::Vector{<:QSym}) where {T<:Number}
        if isequal(arg_c, 0.0)
            return new{T}(0.0, QSym[])
        else
            return new{T}(arg_c, args_nc)
        end
    end
    QMul(args_nc::Vector{<:QSym}) = new{Float64}(1.0, args_nc)
    QMul(s) = new{Float64}(1.0, [s])
    QMul{T}(s) where {T} = new{T}(T(1.0), [s])
    QMul() = new{Float64}(0.0, QSym[])
    QMul{T}() where {T} = new{T}(T(0.0), QSym[])
end
# QMul(x::T) where {T} = @error "QMul constructor does not accept $T"
Base.promote_rule(::Type{QMul{S}}, ::Type{QMul{T}}) where {S,T} = QMul{promote_rule(S, T)}
function Base.convert(::Type{QMul{T}}, x::QMul{S}) where {T<:Number,S<:Number}
    return QMul(convert(T, x.arg_c), x.args_nc)
end

SymbolicUtils.operation(::QMul) = (*)
"""
    arguments(a::QMul)

Return the vector of the factors of [`QMul`](@ref).
"""
SymbolicUtils.arguments(a::QMul) = vcat(a.arg_c, a.args_nc)

function SymbolicUtils.maketerm(::Type{<:QMul}, ::typeof(*), args, metadata)
    args_c = filter(x -> !(x isa QField), args)
    args_nc = filter(x -> x isa QField, args)
    arg_c = isempty(args_c) ? 1 : *(args_c...)
    # isempty(args_nc) && return arg_c
    return QMul(arg_c, args_nc)
end

SymbolicUtils.metadata(a::QMul) = nothing

function Base.adjoint(q::QMul)
    args_nc = map(adjoint, q.args_nc)
    reverse!(args_nc) # TODO fields switch under adjoint right?
    sort!(args_nc; by=position)
    sort!(args_nc; by=ladder)
    return QMul(conj(q.arg_c), args_nc)
end

function Base.isequal(a::QMul, b::QMul)
    isequal(a.arg_c, b.arg_c) || return false
    length(a.args_nc) == length(b.args_nc) || return false
    for (arg_a, arg_b) in zip(a.args_nc, b.args_nc)
        isequal(arg_a, arg_b) || return false
    end
    return true
end

Base.zero(q::QMul) = QMul(0.0, QSym[])
Base.iszero(q::QMul) = iszero(q.arg_c)

"""
    isbulk(q::QTerm)

Checks if a term is in the bulk. A term is bulk if it has no `In` or `Out` position fields ([`AbstractPosition`](@ref)).
"""
isbulk(q::QMul) = all(isbulk.(q.args_nc))
allfields(q::QMul) = q.args_nc

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
        new{Float64}([QMul(s) for s in args])
    end
    QAdd{T}() where {T} = new{T}([QMul{T}()])
end
# QAdd(x::T) where {T} = @error "QAdd constructor does not accept $T"

SymbolicUtils.operation(::QAdd) = (+)
"""
    arguments(a::QAdd)

Return the vector of the arguments of [`QAdd`](@ref).
"""
SymbolicUtils.arguments(a::QAdd) = a.arguments
SymbolicUtils.maketerm(::Type{<:QAdd}, ::typeof(+), args, metadata) = QAdd(args)
SymbolicUtils.metadata(a::QAdd) = nothing

function Base.isequal(a::QAdd, b::QAdd)
    length(arguments(a)) == length(arguments(b)) || return false
    for (arg_a, arg_b) in zip(arguments(a), arguments(b))
        isequal(arg_a, arg_b) || return false
    end
    return true
end

Base.adjoint(q::QAdd) = QAdd(map(adjoint, arguments(q)))
isbulk(q::QAdd) = all(isbulk.(arguments(q)))
allfields(q::QAdd) = reduce(vcat, allfields.(SymbolicUtils.arguments(q)))

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
