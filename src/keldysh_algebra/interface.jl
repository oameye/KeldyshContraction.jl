########################
#   Field Types
########################

"""Abstract type representing expressions involving Keldysh fields."""
abstract type QField end

const CSym = SymbolicUtils.BasicSymbolic{
    isdefined(SymbolicUtils, :SymReal) ? SymbolicUtils.SymReal : Number
}
const SNuN = Union{CSym,Number}
const QSymbol = Union{QField,SNuN}

"""Abstract type representing a fundamental Keldysh field."""
abstract type QSym <: QField end

"""Abstract type representing a noncommutative field expression."""
abstract type QTerm <: QField end

"""Statistics carried statically by every fundamental field."""
abstract type Statistics end

"""Bosonic field statistics."""
struct Boson <: Statistics end

"""Orientation of a path-integral field."""
@enumx Orientation::UInt8 begin
    Unbarred = 0
    Barred = 1
end

"""
Neutral two-valued Keldysh index stored by `Field`.

Bosonic code exposes the semantic aliases `Quantum` and `Classical`. Fermionic
`One`/`Two` aliases can be added later without changing the field representation.
"""
@enumx KeldyshIndex::UInt8 begin
    First = 0
    Second = 1
end

const Quantum = KeldyshIndex.First
const Classical = KeldyshIndex.Second

"""Kind of concrete internal field index."""
@enumx IndexKind::UInt8 begin
    Spin = 0
    Flavor = 1
    Band = 2
    Species = 3
end

"""
A compact concrete internal field index such as spin, flavor, band, or species.

The `Symbol` constructor is provided as a convenient value-level API while the stored
representation uses the one-byte [`IndexKind`](@ref) enum.
"""
struct FieldIndex
    kind::IndexKind.T
    value::Int16
end

function FieldIndex(kind::Symbol, value::Integer)
    k = if kind === :spin
        IndexKind.Spin
    elseif kind === :flavor
        IndexKind.Flavor
    elseif kind === :band
        IndexKind.Band
    elseif kind === :species
        IndexKind.Species
    else
        throw(ArgumentError("unknown field-index kind: $kind"))
    end
    return FieldIndex(k, convert(Int16, value))
end

function Base.isequal(a::FieldIndex, b::FieldIndex)
    return a.kind === b.kind && a.value == b.value
end
Base.hash(x::FieldIndex, h::UInt) = hash(x.kind, hash(x.value, h))
function Base.isless(a::FieldIndex, b::FieldIndex)
    a.kind === b.kind || return Integer(a.kind) < Integer(b.kind)
    return a.value < b.value
end

const MAX_FIELD_INDICES = 4
const EMPTY_FIELD_INDEX = FieldIndex(IndexKind.Spin, Int16(0))

"""
Immutable fixed-capacity storage for field indices.

The number and values of indices are runtime data, while the Julia representation is
fixed and hash-stable. Each slot is compact and isbits; four slots cover the current
spin/flavor/band/species use cases without introducing field-type proliferation.
"""
struct FieldIndices
    data::NTuple{MAX_FIELD_INDICES,FieldIndex}
    n::UInt8
end

FieldIndices() = FieldIndices(ntuple(_ -> EMPTY_FIELD_INDEX, MAX_FIELD_INDICES), UInt8(0))
function FieldIndices(indices::Vararg{FieldIndex,N}) where {N}
    N <= MAX_FIELD_INDICES ||
        throw(ArgumentError("at most $MAX_FIELD_INDICES field indices are supported"))
    data = ntuple(i -> i <= N ? indices[i] : EMPTY_FIELD_INDEX, MAX_FIELD_INDICES)
    return FieldIndices(data, UInt8(N))
end

Base.length(indices::FieldIndices) = Int(indices.n)
function Base.getindex(indices::FieldIndices, i::Int)
    1 <= i <= length(indices) || throw(BoundsError(indices, i))
    return indices.data[i]
end
function Base.iterate(indices::FieldIndices, state::Int=1)
    state > length(indices) && return nothing
    return indices.data[state], state + 1
end
function Base.isequal(a::FieldIndices, b::FieldIndices)
    length(a) == length(b) || return false
    @inbounds for i in 1:length(a)
        isequal(a.data[i], b.data[i]) || return false
    end
    return true
end
function Base.hash(indices::FieldIndices, h::UInt)
    h = hash(FieldIndices, hash(indices.n, h))
    @inbounds for i in 1:length(indices)
        h = hash(indices.data[i], h)
    end
    return h
end
function Base.isless(a::FieldIndices, b::FieldIndices)
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        isequal(a.data[i], b.data[i]) && continue
        return isless(a.data[i], b.data[i])
    end
    return length(a) < length(b)
end

################################
# Interface for SymbolicUtils
################################

TermInterface.head(::QField) = :call
SymbolicUtils.iscall(::QSym) = false
SymbolicUtils.iscall(::QTerm) = true
SymbolicUtils.iscall(::Type{T}) where {T<:QTerm} = true
TermInterface.metadata(::QSym) = nothing

for f in SymbolicUtils.basic_diadic
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)), Ts::Type{<:QField}...) =
        promote_type(Ts...)
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)), T::Type{<:QField}, Ts...) = T
    @eval SymbolicUtils.promote_symtype(
        ::$(typeof(f)), T::Type{<:QField}, S::Type{<:Number}
    ) = T
    @eval SymbolicUtils.promote_symtype(
        ::$(typeof(f)), T::Type{<:Number}, S::Type{<:QField}
    ) = S
    @eval SymbolicUtils.promote_symtype(
        ::$(typeof(f)), T::Type{<:QField}, S::Type{<:QField}
    ) = promote_type(T, S)
end

SymbolicUtils.symtype(x::T) where {T<:QField} = T

Base.isone(::QField) = false
Base.iszero(::QField) = false

#########################
# Field enum properties
#########################

name(ϕ::QSym) = ϕ.name

"""Regularisation enum used for equal-time Keldysh contractions."""
@enumx Regularisation::Int8 begin
    Plus = 1
    Zero = 0
    Minus = -1
end

subtraction(x::NTuple{2,Regularisation.T}) = -(Int.(x)...)
function subtraction(x::Vector{Regularisation.T})
    length(x) == 2 || throw(ArgumentError("regularisation subtraction expects two values"))
    return subtraction(Tuple(x))
end

#########################
#       Position
#########################

"""Concrete position of a Keldysh field."""
struct Position
    index::Int8
end

function Bulk(i::Integer=1)
    i > 0 || throw(ArgumentError("Bulk index must be positive"))
    return Position(convert(Int8, i))
end
In() = Position(typemax(Int8))
Out() = Position(typemin(Int8))

Base.isless(x::Position, y::Position) = x.index < y.index

is_bulk(x::Position) = typemin(Int8) < index(x) < typemax(Int8)
is_in(x::Position) = isequal(index(x), typemax(Int8))
is_out(x::Position) = isequal(index(x), typemin(Int8))

"""Exchange `In()` and `Out()` while leaving bulk coordinates unchanged."""
function swap_in_out(x::Position)
    is_in(x) && return Out()
    is_out(x) && return In()
    return x
end

index(p::Position) = p.index
Base.Int(p::Position) = index(p)

has_in(ps)::Bool = In() ∈ ps
has_out(ps)::Bool = Out() ∈ ps
