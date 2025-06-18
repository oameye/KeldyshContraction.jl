# These has function are defined such that comparison functions work, e.g., `unique`.
Base.hash(q::QMul, h::UInt) = hash(QMul, hash(q.arg_c, SymbolicUtils.hashvec(q.args_nc, h)))

Base.hash(q::QAdd, h::UInt) = hash(QAdd, SymbolicUtils.hashvec(arguments(q), h))

function Base.hash(h::Union{KeldyshContour.T,Regularisation.T}, i::UInt)
    return hash(Int(h), i)
end

for f in [:Destroy, :Create]
    @eval function Base.hash(op::T, h::UInt) where {T<:($(f))}
        return hash(
            $(f),
            hash(
                name(op), hash(contour(op), hash(position(op), hash(regularisation(op), h)))
            ),
        )
    end
end
