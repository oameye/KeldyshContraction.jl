struct Momentum
    index::Int8
    function Momentum(index)
        index < 0 && error("Momentum index must be non-negative, got $index")
        return new(index)
    end
end
Base.isequal(m1::Momentum, m2::Momentum) = m1.index == m2.index
struct Momenta
    prefactors::Vector{Int}
    Momenta::Vector{Momentum}
    function Momenta(prefactors::Vector{Int}, Momenta::Vector{Momentum})
        @assert length(prefactors) == length(Momenta) "Length of prefactors and momenta must match"
        idxs = findall(x -> !iszero(x), prefactors)
        return new(prefactors[idxs], Momenta[idxs])
    end
end
function Momenta(idx::Int)
    return Momenta([1], [Momentum(idx)])
end
function Base.isequal(m1::Momenta, m2::Momenta)
    isequal(m1.prefactors, m2.prefactors) && isequal(m1.Momenta, m2.Momenta)
end
# SmallCollections.default(::Type{Momenta}) = Momenta(Int[], Momentum[])
