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
Index-kind names in enum order, so `INDEX_KIND_NAMES[Int(kind) + 1]` names a kind and
`findfirst` maps a name back. One table replaces the per-kind branches everywhere below.
"""
const INDEX_KIND_NAMES = (:spin, :flavor, :band, :species)

const MAX_FIELD_INDICES = length(INDEX_KIND_NAMES)
const NO_FIELD_INDEX = typemin(Int16)

"""
$(DocStringExtensions.TYPEDEF)

A compact concrete internal field index such as spin, flavor, band, or species.

The `Symbol` constructor is a convenient value-level API while the stored representation
uses a one-byte index-kind enum.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FieldIndex
    "Which kind of internal index this is"
    kind::IndexKind.T
    "The index value"
    value::Int16
end

function FieldIndex(kind::Symbol, value::Integer)
    slot = findfirst(==(kind), INDEX_KIND_NAMES)
    slot === nothing && throw(ArgumentError("unknown field-index kind: $kind"))
    return FieldIndex(IndexKind.T(slot - 1), convert(Int16, value))
end

"""Position of a kind in the [`FieldIndices`](@ref) slot tuple."""
@inline _slot(kind::IndexKind.T) = Int(kind) + 1

Base.isequal(a::FieldIndex, b::FieldIndex) = a.kind === b.kind && a.value == b.value
Base.hash(x::FieldIndex, h::UInt) = hash(x.kind, hash(x.value, h))
function Base.isless(a::FieldIndex, b::FieldIndex)
    a.kind === b.kind || return Int(a.kind) < Int(b.kind)
    return a.value < b.value
end

"""
$(DocStringExtensions.TYPEDEF)

Compact value-level storage for the supported field-index kinds.

Each kind has one `Int16` slot; `typemin(Int16)` denotes absence. This keeps every
`Field{S}` concrete and compact without making index values or index presence type
parameters. Iteration yields present [`FieldIndex`](@ref) values in canonical
spin/flavor/band/species order.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FieldIndices
    "Spin index, or `NO_FIELD_INDEX` when absent"
    spin::Int16
    "Flavor index, or `NO_FIELD_INDEX` when absent"
    flavor::Int16
    "Band index, or `NO_FIELD_INDEX` when absent"
    band::Int16
    "Species index, or `NO_FIELD_INDEX` when absent"
    species::Int16
end

"""Slot values in canonical kind order, so every query below is one tuple lookup."""
@inline slots(x::FieldIndices) = (x.spin, x.flavor, x.band, x.species)

@inline _has_index(value::Int16) = value != NO_FIELD_INDEX

FieldIndices() = FieldIndices(ntuple(_ -> NO_FIELD_INDEX, MAX_FIELD_INDICES)...)

function FieldIndices(indices::Vararg{FieldIndex,N}) where {N}
    N <= MAX_FIELD_INDICES ||
        throw(ArgumentError("at most $MAX_FIELD_INDICES field indices are supported"))

    values = ntuple(_ -> NO_FIELD_INDEX, MAX_FIELD_INDICES)
    for idx in indices
        idx.value == NO_FIELD_INDEX && throw(
            ArgumentError("$(NO_FIELD_INDEX) is reserved as the absent-index sentinel")
        )
        slot = _slot(idx.kind)
        _has_index(values[slot]) &&
            throw(ArgumentError("duplicate $(INDEX_KIND_NAMES[slot]) index"))
        values = Base.setindex(values, idx.value, slot)
    end
    return FieldIndices(values...)
end

Base.length(indices::FieldIndices) = count(_has_index, slots(indices))

function Base.getindex(indices::FieldIndices, i::Int)
    1 <= i <= length(indices) || throw(BoundsError(indices, i))
    n = 0
    for (slot, value) in enumerate(slots(indices))
        _has_index(value) || continue
        n += 1
        n == i && return FieldIndex(IndexKind.T(slot - 1), value)
    end
    return throw(BoundsError(indices, i))
end
function Base.iterate(indices::FieldIndices, state::Int=1)
    state > length(indices) && return nothing
    return indices[state], state + 1
end
Base.eltype(::Type{FieldIndices}) = FieldIndex
Base.IteratorSize(::Type{FieldIndices}) = Base.HasLength()

Base.isequal(a::FieldIndices, b::FieldIndices) = slots(a) == slots(b)
Base.hash(indices::FieldIndices, h::UInt) = hash(slots(indices), h)
Base.isless(a::FieldIndices, b::FieldIndices) = slots(a) < slots(b)

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

"""
$(DocStringExtensions.TYPEDEF)

Concrete position of a Keldysh field. A position has three cases:
- [`In`](@ref): index is `typemax(Int8)`, representing the incoming external field;
- [`Out`](@ref): index is `typemin(Int8)`, representing the outgoing external field;
- [`Bulk`](@ref): index is a positive integer naming an interaction vertex.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Position
    "Encoded coordinate; the sentinel extremes mark the external legs"
    index::Int8
end

"""
$(DocStringExtensions.SIGNATURES)

Create a `Bulk` position with positive integer index `i`.
"""
function Bulk(i::Integer=1)
    i > 0 || throw(ArgumentError("Bulk index must be positive"))
    return Position(convert(Int8, i))
end

"""
$(DocStringExtensions.SIGNATURES)

Create a `Position` representing the incoming external field.
"""
In() = Position(typemax(Int8))

"""
$(DocStringExtensions.SIGNATURES)

Create a `Position` representing the outgoing external field.
"""
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
