"""One exact statistical distribution factor `F_family(momentum)` carried through frequency reduction."""
struct StatisticalAtom{S<:Statistics}
    family::FieldFamily{S}
    momentum::LinearMomentum
end

statistics(::StatisticalAtom{S}) where {S<:Statistics} = S
statistical_family(atom::StatisticalAtom) = atom.family
momentum(atom::StatisticalAtom) = atom.momentum

function Base.isequal(a::StatisticalAtom{S}, b::StatisticalAtom{S}) where {S<:Statistics}
    return isequal(a.family, b.family) && isequal(a.momentum, b.momentum)
end
Base.:(==)(a::StatisticalAtom, b::StatisticalAtom) = isequal(a, b)
function Base.hash(atom::StatisticalAtom, h::UInt)
    return hash(atom.momentum, hash(atom.family, hash(StatisticalAtom, h)))
end
function Base.isless(a::StatisticalAtom{S}, b::StatisticalAtom{S}) where {S<:Statistics}
    isequal(a.family, b.family) || return isless(a.family, b.family)
    isequal(a.momentum, b.momentum) && return false
    return _linear_momentum_isless(a.momentum, b.momentum)
end

"""Canonical commutative product of exact statistical distribution factors."""
struct StatisticalMonomial{S<:Statistics}
    factors::Vector{StatisticalAtom{S}}

    function StatisticalMonomial{S}(
        factors::Vector{StatisticalAtom{S}}, ::Val{:canonical}
    ) where {S<:Statistics}
        return new{S}(factors)
    end
end

function StatisticalMonomial{S}() where {S<:Statistics}
    return StatisticalMonomial{S}(StatisticalAtom{S}[], Val(:canonical))
end

function StatisticalMonomial(
    factors::AbstractVector{StatisticalAtom{S}}
) where {S<:Statistics}
    canonical = collect(factors)
    sort!(canonical)
    return StatisticalMonomial{S}(canonical, Val(:canonical))
end

Base.length(monomial::StatisticalMonomial) = length(monomial.factors)
Base.isempty(monomial::StatisticalMonomial) = isempty(monomial.factors)
Base.iterate(monomial::StatisticalMonomial) = iterate(monomial.factors)
Base.iterate(monomial::StatisticalMonomial, state) = iterate(monomial.factors, state)
Base.eltype(::Type{StatisticalMonomial{S}}) where {S} = StatisticalAtom{S}
Base.IteratorSize(::Type{<:StatisticalMonomial}) = Base.HasLength()

function Base.isequal(a::StatisticalMonomial{S}, b::StatisticalMonomial{S}) where {S}
    return isequal(a.factors, b.factors)
end
Base.:(==)(a::StatisticalMonomial, b::StatisticalMonomial) = isequal(a, b)
function Base.hash(monomial::StatisticalMonomial, h::UInt)
    return hash(monomial.factors, hash(StatisticalMonomial, h))
end
function Base.isless(a::StatisticalMonomial{S}, b::StatisticalMonomial{S}) where {S}
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        ai = a.factors[i]
        bi = b.factors[i]
        isequal(ai, bi) || return isless(ai, bi)
    end
    return length(a) < length(b)
end

function Base.:*(a::StatisticalMonomial{S}, b::StatisticalMonomial{S}) where {S}
    return StatisticalMonomial(vcat(a.factors, b.factors))
end
