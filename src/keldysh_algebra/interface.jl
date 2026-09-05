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

The `Symbol` constructor is a convenient value-level API while the stored representation
uses a one-byte index-kind enum.
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
    a.kind === b.kind || return Int(a.kind) < Int(b.kind)
    return a.value < b.value
end

const MAX_FIELD_INDICES = 4
const NO_FIELD_INDEX = typemin(Int16)

"""
Compact value-level storage for the supported field-index kinds.

Each kind has one `Int16` slot; `typemin(Int16)` denotes absence. This keeps every
`Field{S}` concrete and compact without making index values or index presence type
parameters. Iteration yields present [`FieldIndex`](@ref) values in canonical
spin/flavor/band/species order.
"""
struct FieldIndices
    spin::Int16
    flavor::Int16
    band::Int16
    species::Int16
end

FieldIndices() = FieldIndices(NO_FIELD_INDEX, NO_FIELD_INDEX, NO_FIELD_INDEX, NO_FIELD_INDEX)
function FieldIndices(indices::Vararg{FieldIndex,N}) where {N}
    N <= MAX_FIELD_INDICES ||
        throw(ArgumentError("at most $MAX_FIELD_INDICES field indices are supported"))

    spin = NO_FIELD_INDEX
    flavor = NO_FIELD_INDEX
    band = NO_FIELD_INDEX
    species = NO_FIELD_INDEX

    for idx in indices
        idx.value == NO_FIELD_INDEX &&
            throw(ArgumentError("$(NO_FIELD_INDEX) is reserved as the absent-index sentinel"))
        if idx.kind === IndexKind.Spin
            spin == NO_FIELD_INDEX || throw(ArgumentError("duplicate spin index"))
            spin = idx.value
        elseif idx.kind === IndexKind.Flavor
            flavor == NO_FIELD_INDEX || throw(ArgumentError("duplicate flavor index"))
            flavor = idx.value
        elseif idx.kind === IndexKind.Band
            band == NO_FIELD_INDEX || throw(ArgumentError("duplicate band index"))
            band = idx.value
        else
            species == NO_FIELD_INDEX || throw(ArgumentError("duplicate species index"))
            species = idx.value
        end
    end

    return FieldIndices(spin, flavor, band, species)
end

@inline _has_index(value::Int16) = value != NO_FIELD_INDEX
function Base.length(indices::FieldIndices)
    return Int(_has_index(indices.spin)) +
           Int(_has_index(indices.flavor)) +
           Int(_has_index(indices.band)) +
           Int(_has_index(indices.species))
end

function Base.getindex(indices::FieldIndices, i::Int)
    1 <= i <= length(indices) || throw(BoundsError(indices, i))
    n = 0
    if _has_index(indices.spin)
        n += 1
        n == i && return FieldIndex(IndexKind.Spin, indices.spin)
    end
    if _has_index(indices.flavor)
        n += 1
        n == i && return FieldIndex(IndexKind.Flavor, indices.flavor)
    end
    if _has_index(indices.band)
        n += 1
        n == i && return FieldIndex(IndexKind.Band, indices.band)
    end
    return FieldIndex(IndexKind.Species, indices.species)
end
function Base.iterate(indices::FieldIndices, state::Int=1)
    state > length(indices) && return nothing
    return indices[state], state + 1
end

function Base.isequal(a::FieldIndices, b::FieldIndices)
    return a.spin == b.spin &&
           a.flavor == b.flavor &&
           a.band == b.band &&
           a.species == b.species
end
Base.hash(indices::FieldIndices, h::UInt) =
    hash(indices.species, hash(indices.band, hash(indices.flavor, hash(indices.spin, h))))
function Base.isless(a::FieldIndices, b::FieldIndices)
    return (a.spin, a.flavor, a.band, a.species) < (b.spin, b.flavor, b.band, b.species)
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
