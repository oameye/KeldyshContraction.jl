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

"""Statistics carried by a Keldysh field."""
abstract type Statistics end

"""Bosonic field statistics."""
struct Boson <: Statistics end

"""Orientation of a path-integral field."""
@enumx Orientation::UInt8 begin
    Unbarred = 0
    Barred = 1
end

"""Two-valued Keldysh index stored by `Field`."""
@enumx KeldyshIndex::UInt8 begin
    First = 0
    Second = 1
end

const Quantum = KeldyshIndex.First
const Classical = KeldyshIndex.Second

"""Kind of internal field index."""
@enumx IndexKind::UInt8 begin
    Spin = 0
    Flavor = 1
    Band = 2
    Species = 3
end

const INDEX_KIND_NAMES = (:spin, :flavor, :band, :species)
const MAX_FIELD_INDICES = length(INDEX_KIND_NAMES)
const NO_FIELD_INDEX = typemin(Int16)

"""Internal field index such as spin, flavor, band, or species."""
struct FieldIndex
    kind::IndexKind.T
    value::Int16
end

function FieldIndex(kind::Symbol, value::Integer)
    slot = findfirst(==(kind), INDEX_KIND_NAMES)
    slot === nothing && throw(ArgumentError("unknown field-index kind: $kind"))
    return FieldIndex(IndexKind.T(slot - 1), convert(Int16, value))
end

@inline field_index_slot(kind::IndexKind.T) = Int(kind) + 1

Base.isequal(a::FieldIndex, b::FieldIndex) = a.kind === b.kind && a.value == b.value
Base.hash(x::FieldIndex, h::UInt) = hash(x.kind, hash(x.value, h))
function Base.isless(a::FieldIndex, b::FieldIndex)
    a.kind === b.kind || return Int(a.kind) < Int(b.kind)
    return a.value < b.value
end

"""Internal indices of a Keldysh field."""
struct FieldIndices
    spin::Int16
    flavor::Int16
    band::Int16
    species::Int16
end

@inline slots(x::FieldIndices) = (x.spin, x.flavor, x.band, x.species)
@inline has_field_index(value::Int16) = value != NO_FIELD_INDEX

FieldIndices() = FieldIndices(ntuple(_ -> NO_FIELD_INDEX, MAX_FIELD_INDICES)...)

function FieldIndices(indices::Vararg{FieldIndex,N}) where {N}
    N <= MAX_FIELD_INDICES ||
        throw(ArgumentError("at most $MAX_FIELD_INDICES field indices are supported"))

    values = ntuple(_ -> NO_FIELD_INDEX, MAX_FIELD_INDICES)
    for idx in indices
        idx.value == NO_FIELD_INDEX && throw(
            ArgumentError("$(NO_FIELD_INDEX) is reserved as the absent-index sentinel")
        )
        slot = field_index_slot(idx.kind)
        has_field_index(values[slot]) &&
            throw(ArgumentError("duplicate $(INDEX_KIND_NAMES[slot]) index"))
        values = Base.setindex(values, idx.value, slot)
    end
    return FieldIndices(values...)
end

Base.length(indices::FieldIndices) = count(has_field_index, slots(indices))

function Base.getindex(indices::FieldIndices, i::Int)
    1 <= i <= length(indices) || throw(BoundsError(indices, i))
    n = 0
    for (slot, value) in enumerate(slots(indices))
        has_field_index(value) || continue
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

"""Regularisation used for equal-time Keldysh contractions."""
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

"""Position of a Keldysh field."""
struct Position
    index::Int8
end

"""Create a bulk position with index `i`."""
function Bulk(i::Integer=1)
    i > 0 || throw(ArgumentError("Bulk index must be positive"))
    return Position(convert(Int8, i))
end

"""Incoming external position."""
In() = Position(typemax(Int8))

"""Outgoing external position."""
Out() = Position(typemin(Int8))

Base.isless(x::Position, y::Position) = x.index < y.index

is_bulk(x::Position) = typemin(Int8) < index(x) < typemax(Int8)
is_in(x::Position) = isequal(index(x), typemax(Int8))
is_out(x::Position) = isequal(index(x), typemin(Int8))

"""Exchange `In()` and `Out()` while leaving bulk positions unchanged."""
function swap_in_out(x::Position)
    is_in(x) && return Out()
    is_out(x) && return In()
    return x
end

index(p::Position) = p.index
Base.Int(p::Position) = index(p)

has_in(ps)::Bool = In() ∈ ps
has_out(ps)::Bool = Out() ∈ ps
