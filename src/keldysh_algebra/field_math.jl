# *,+,/,- math of the QField and QSym types

function Base.:*(a::QSym, b::QSym)
    args = QSym[a, b]
    sort!(args; by=position)
    sort!(args; by=ladder)
    return QMul(1, args)
end

function Base.:*(a::QSym, b::Number)
    # SymbolicUtils._iszero(b) && return b
    # SymbolicUtils._isone(b) && return a
    return QMul(b, QSym[a])
end
Base.:*(b::Number, a::QField) = a * b

function Base.:*(a::QMul, b::Number)
    # SymbolicUtils._iszero(b) && return b
    # SymbolicUtils._isone(b) && return a
    arg_c = a.arg_c * b
    return QMul(arg_c, a.args_nc)
end

function Base.:*(a::QSym, b::QMul)
    args_nc = vcat(a, b.args_nc)
    sort!(args_nc; by=position)
    sort!(args_nc; by=ladder)
    return QMul(b.arg_c, args_nc)
end
Base.:*(a::QMul, b::QSym) = b * a

function Base.:*(a::QMul, b::QMul)
    args_nc = vcat(a.args_nc, b.args_nc)
    sort!(args_nc; by=position)
    sort!(args_nc; by=ladder)
    arg_c = a.arg_c * b.arg_c
    return QMul(arg_c, args_nc)
end

Base.:/(a::QField, b::Number) = (1 // b) * a ## TODO: maybe rationalize?

## Powers
function Base.:^(a::QField, n::Integer)
    # Type-stable implementation using Val
    iszero(n) && return QMul(0, QSym[])
    result = QMul(1, QSym[a])
    for _ in 2:n
        result *= a
    end
    return result
end

## Addition

Base.:-(a::QField) = -1 * a
Base.:-(a::Number, b::QField) = a + (-b)
Base.:-(a::QField, b::Number) = a + (-b)
Base.:-(a::QField, b::QField) = a + (-b)

function Base.:+(a::QField, b::Number)
    SymbolicUtils._iszero(b) && return QAdd([a])
    return QAdd([a, b])
end
Base.:+(a::Number, b::QField) = +(b, a)
function Base.:+(a::QAdd, b::Number)
    SymbolicUtils._iszero(b) && return a
    args = vcat(arguments(a), b)
    return QAdd(args)
end

function Base.:+(a::QSym, b::QSym)
    args = [a, b]
    return QAdd(args)
end
function Base.:+(a::QMul{T}, b::QMul{S}) where {T,S}
    TT = promote_type(T, S)
    args = QMul{TT}[a, b]
    return QAdd(args)
end
function Base.:+(a::QMul{T}, b::QSym) where {T}
    args = QMul{T}[a, QMul{T}(b)]
    return QAdd(args)
end
function Base.:+(b::QSym, a::QMul{T}) where {T}
    args = QMul{T}[a, QMul{T}(b)]
    return QAdd(args)
end

function Base.:+(a::QMul{T}, b::QAdd{S}) where {T,S}
    TT = promote_type(T, S)
    args = QMul{TT}[_b for _b in arguments(b)]
    return QAdd(push!(args, a))
end
function Base.:+(b::QAdd{S}, a::QMul{T}) where {T,S}
    TT = promote_type(T, S)
    args = QMul{TT}[_b for _b in arguments(b)]
    return QAdd(push!(args, a))
end

function Base.:+(a::QAdd{T}, b::QSym) where {T}
    args = arguments(a)
    l = length(args)
    args′ = QMul{T}[i>l ? QMul{T}(b) : args[i] for i in 1:(l + 1)]
    return QAdd(args′)
end
function Base.:+(b::QSym, a::QAdd{T}) where {T}
    args = arguments(a)
    l = length(args)
    args′ = QMul{T}[i>l ? QMul{T}(b) : args[i] for i in 1:(l + 1)]
    return QAdd(args′)
end
function Base.:+(a::QAdd{T}, b::QAdd{S}) where {T,S}
    TT = promote_type(T, S)
    args1 = arguments(a)
    args2 = arguments(b)
    l1 = length(args1)
    l2 = length(args2)
    args′ = QMul{TT}[i>l1 ? args2[i - l1] : args1[i] for i in 1:(l1 + l2)]
    return QAdd(args′)
end

function Base.:*(a::QAdd{T}, b::S) where {T,S<:Number}
    TT = promote_type(T, S)
    args = QMul{TT}[a_ * b for a_ in arguments(a)]
    flatten_adds!(args)
    isempty(args) && return QAdd{TT}()
    return QAdd(args)
end

function Base.:*(a::QMul{S}, b::QAdd{T}) where {T,S}
    TT = promote_type(T, S)
    args = QMul{TT}[b_ * a for b_ in arguments(b)]
    flatten_adds!(args)
    isempty(args) && return QAdd{TT}()
    return QAdd(args)
end
function Base.:*(b::QAdd{T}, a::QMul{S}) where {T,S}
    TT = promote_type(T, S)
    args = QMul{TT}[b_ * a for b_ in arguments(b)]
    flatten_adds!(args)
    isempty(args) && return QAdd{TT}()
    return QAdd(args)
end

function Base.:*(a::QSym, b::QAdd{T}) where {T}
    args = QMul{T}[b_ * a for b_ in arguments(b)]
    flatten_adds!(args)
    isempty(args) && return QAdd{T}()
    return QAdd(args)
end
function Base.:*(b::QAdd{T}, a::QSym) where {T}
    args = QMul{T}[b_ * a for b_ in arguments(b)]
    flatten_adds!(args)
    isempty(args) && return QAdd{T}()
    return QAdd(args)
end

function Base.:*(a::QAdd{S}, b::QAdd{T}) where {T,S}
    TT = promote_type(T, S)
    args = QMul{TT}[]
    for a_ in arguments(a), b_ in arguments(b)
        push!(args, a_ * b_)
    end
    flatten_adds!(args)
    isempty(args) && return QAdd{TT}()
    return QAdd(args)
end

function flatten_adds!(args)
    i = 1
    while i <= length(args)
        x = args[i]
        if x isa QAdd
            append!(args, arguments(x))
            deleteat!(args, i)
        elseif SymbolicUtils._iszero(x) || isequal(x, 0)
            deleteat!(args, i)
        else
            i += 1
        end
    end
    return args
end
