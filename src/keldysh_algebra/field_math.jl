# Arithmetic for Field{S}, QMul{C,S}, and QAdd{C,S}.

@inline function _converted_mul(::Type{C}, q::QMul{D,S}) where {C<:Number,D<:Number,S}
    return convert(QMul{C,S}, q)
end

function _normalized_add(args::Vector{QMul{C,S}}) where {C<:Number,S<:Statistics}
    out = QMul{C,S}[term for term in args if !iszero(term)]
    isempty(out) && push!(out, QMul{C,S}(zero(C), Field{S}[]))
    return QAdd{C,S}(out)
end

########################
# Multiplication
########################

function Base.:*(a::Field{S}, b::Field{S}) where {S<:Statistics}
    return QMul(1, Field{S}[a, b])
end

Base.:*(a::Field{S}, b::C) where {S<:Statistics,C<:Number} = QMul(b, Field{S}[a])
Base.:*(b::Number, a::Field) = a * b

function Base.:*(a::QMul{C,S}, b::D) where {C<:Number,D<:Number,S<:Statistics}
    coeff = a.arg_c * b
    return QMul(coeff, copy(a.args_nc))
end
Base.:*(b::Number, a::QMul) = a * b

function Base.:*(a::Field{S}, b::QMul{C,S}) where {C<:Number,S<:Statistics}
    return QMul(b.arg_c, vcat(Field{S}[a], b.args_nc))
end
Base.:*(a::QMul{C,S}, b::Field{S}) where {C<:Number,S<:Statistics} = b * a

function Base.:*(a::QMul{C1,S}, b::QMul{C2,S}) where {C1,C2,S<:Statistics}
    coeff = a.arg_c * b.arg_c
    return QMul(coeff, vcat(a.args_nc, b.args_nc))
end

Base.:/(a::Field, b::Number) = a * inv(b)
Base.:/(a::QMul, b::Number) = a * inv(b)
Base.:/(a::QAdd, b::Number) = a * inv(b)
Base.://(a::Field, b::Number) = (1 // b) * a
Base.://(a::QMul, b::Number) = (1 // b) * a
Base.://(a::QAdd, b::Number) = (1 // b) * a

########################
# Powers
########################

function Base.:^(a::Field{S}, n::Integer) where {S<:Statistics}
    n < 0 && throw(DomainError(n, "negative field powers are not supported"))
    result = QMul(1, Field{S}[])
    for _ in 1:n
        result *= a
    end
    return result
end

function Base.:^(a::QMul{C,S}, n::Integer) where {C<:Number,S<:Statistics}
    n < 0 && throw(DomainError(n, "negative field powers are not supported"))
    result = QMul{C,S}(one(C), Field{S}[])
    for _ in 1:n
        result *= a
    end
    return result
end

function Base.:^(a::QAdd{C,S}, n::Integer) where {C<:Number,S<:Statistics}
    n < 0 && throw(DomainError(n, "negative field powers are not supported"))
    result = one(a)
    for _ in 1:n
        result *= a
    end
    return result
end

########################
# Addition / subtraction
########################

Base.:-(a::Field) = -1 * a
Base.:-(a::QMul) = -1 * a
Base.:-(a::QAdd) = -1 * a

Base.:-(a::Number, b::QField) = a + (-b)
Base.:-(a::QField, b::Number) = a + (-b)
Base.:-(a::QField, b::QField) = a + (-b)

function Base.:+(a::Field{S}, b::D) where {S<:Statistics,D<:Number}
    C = promote_type(Int, D)
    terms = QMul{C,S}[QMul(convert(C, 1), Field{S}[a]), QMul(convert(C, b), Field{S}[])]
    return _normalized_add(terms)
end
Base.:+(b::Number, a::Field) = a + b

function Base.:+(a::QMul{C,S}, b::D) where {C<:Number,D<:Number,S<:Statistics}
    P = promote_type(C, D)
    terms = QMul{P,S}[_converted_mul(P, a), QMul(convert(P, b), Field{S}[])]
    return _normalized_add(terms)
end
Base.:+(b::Number, a::QMul) = a + b

function Base.:+(a::QAdd{C,S}, b::D) where {C<:Number,D<:Number,S<:Statistics}
    P = promote_type(C, D)
    args = QMul{P,S}[_converted_mul(P, term) for term in a.arguments]
    push!(args, QMul(convert(P, b), Field{S}[]))
    return _normalized_add(args)
end
Base.:+(b::Number, a::QAdd) = a + b

function Base.:+(a::Field{S}, b::Field{S}) where {S<:Statistics}
    return _normalized_add(QMul{Int,S}[QMul(a), QMul(b)])
end

function Base.:+(a::QMul{C1,S}, b::QMul{C2,S}) where {C1,C2,S<:Statistics}
    P = promote_type(C1, C2)
    return _normalized_add(QMul{P,S}[_converted_mul(P, a), _converted_mul(P, b)])
end

function Base.:+(a::QMul{C,S}, b::Field{S}) where {C<:Number,S<:Statistics}
    return _normalized_add(QMul{C,S}[a, QMul(convert(C, 1), Field{S}[b])])
end
Base.:+(b::Field{S}, a::QMul{C,S}) where {C<:Number,S<:Statistics} = a + b

function Base.:+(a::QAdd{C,S}, b::Field{S}) where {C<:Number,S<:Statistics}
    args = copy(a.arguments)
    push!(args, QMul(convert(C, 1), Field{S}[b]))
    return _normalized_add(args)
end
Base.:+(b::Field{S}, a::QAdd{C,S}) where {C<:Number,S<:Statistics} = a + b

function Base.:+(a::QMul{C1,S}, b::QAdd{C2,S}) where {C1,C2,S<:Statistics}
    P = promote_type(C1, C2)
    args = QMul{P,S}[_converted_mul(P, term) for term in b.arguments]
    push!(args, _converted_mul(P, a))
    return _normalized_add(args)
end
Base.:+(b::QAdd{C2,S}, a::QMul{C1,S}) where {C1,C2,S<:Statistics} = a + b

function Base.:+(a::QAdd{C1,S}, b::QAdd{C2,S}) where {C1,C2,S<:Statistics}
    P = promote_type(C1, C2)
    args = QMul{P,S}[]
    append!(args, (_converted_mul(P, term) for term in a.arguments))
    append!(args, (_converted_mul(P, term) for term in b.arguments))
    return _normalized_add(args)
end

########################
# Distributive products
########################

function Base.:*(a::QAdd{C,S}, b::D) where {C<:Number,D<:Number,S<:Statistics}
    P = Base.promote_op(*, C, D)
    args = QMul{P,S}[term * b for term in a.arguments]
    return _normalized_add(args)
end
Base.:*(b::Number, a::QAdd) = a * b

function Base.:*(a::QMul{C1,S}, b::QAdd{C2,S}) where {C1,C2,S<:Statistics}
    P = Base.promote_op(*, C1, C2)
    return _normalized_add(QMul{P,S}[a * term for term in b.arguments])
end
Base.:*(b::QAdd{C2,S}, a::QMul{C1,S}) where {C1,C2,S<:Statistics} = a * b

function Base.:*(a::Field{S}, b::QAdd{C,S}) where {C<:Number,S<:Statistics}
    return _normalized_add(QMul{C,S}[a * term for term in b.arguments])
end
Base.:*(b::QAdd{C,S}, a::Field{S}) where {C<:Number,S<:Statistics} = a * b

function Base.:*(a::QAdd{C1,S}, b::QAdd{C2,S}) where {C1,C2,S<:Statistics}
    P = Base.promote_op(*, C1, C2)
    args = QMul{P,S}[]
    for ta in a.arguments, tb in b.arguments
        push!(args, ta * tb)
    end
    return _normalized_add(args)
end

"""Remove zero terms while preserving the concrete element type."""
function flatten_adds!(args::Vector{QMul{C,S}}) where {C,S}
    filter!(x -> !iszero(x), args)
    return args
end
