# Hash functions are defined consistently with `isequal` so fields/terms remain safe
# dictionary keys during canonicalization.
_hashvec(xs, h::UInt) = foldr(hash, xs; init=h)

function Base.hash(q::QMul, h::UInt)
    isempty(q.args_nc) && return hash(q.arg_c, h)
    length(q.args_nc) == 1 && isone(q.arg_c) && return hash(first(q.args_nc), h)
    return hash(QMul, hash(q.arg_c, _hashvec(q.args_nc, h)))
end

function Base.hash(q::QAdd, h::UInt)
    if length(q.arguments) == 1
        term = first(q.arguments)
        length(term.args_nc) == 1 && isone(term.arg_c) && return hash(term, h)
    end
    return hash(QAdd, _hashvec(q.arguments, h))
end

function Base.hash(h::Union{KeldyshIndex.T,Orientation.T,Regularisation.T}, i::UInt)
    return hash(Int(h), i)
end

function Base.hash(f::FieldFamily{S}, h::UInt) where {S<:Statistics}
    return hash(FieldFamily{S}, hash(name(f), hash(field_indices(f), h)))
end

function Base.hash(f::Field{S}, h::UInt) where {S<:Statistics}
    return hash(
        Field{S},
        hash(
            field_family(f),
            hash(
                orientation(f),
                hash(
                    keldysh_index(f),
                    hash(position(f), hash(regularisation(f), h)),
                ),
            ),
        ),
    )
end

Base.hash(p::Position, h::UInt) = hash(Position, hash(p.index, h))
