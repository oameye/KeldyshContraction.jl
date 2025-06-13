########################
#   Field Types
########################

"""
    QField

Abstract type representing any expression involving Fields.
"""
abstract type QField end
const CSym = SymbolicUtils.Symbolic{<:Number}
const SNuN = Union{CSym,Number}
const QSymbol = Union{QField,SNuN}

"""
    QSym <: QField

Abstract type representing fundamental Field types.
"""
abstract type QSym <: QField end

"""
    QTerm <: QField

Abstract type representing noncommutative expressions.
"""
abstract type QTerm <: QField end

################################
# Interface for SymbolicUtils
################################

TermInterface.head(::QField) = :call
SymbolicUtils.iscall(::QSym) = false
SymbolicUtils.iscall(::QTerm) = true
SymbolicUtils.iscall(::Type{T}) where {T<:QTerm} = true
TermInterface.metadata(x::QSym) = nothing

# Symbolic type promotion
for f in SymbolicUtils.basic_diadic # [+, -, *, /, //, \, ^]
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)), Ts::Type{<:QField}...) = promote_type(
        Ts...
    )
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)), T::Type{<:QField}, Ts...) = T
    @eval SymbolicUtils.promote_symtype(
        ::$(typeof(f)), T::Type{<:QField}, S::Type{<:Number}
    ) = T
    @eval SymbolicUtils.promote_symtype(
        ::$(typeof(f)), T::Type{<:Number}, S::Type{<:QField}
    ) = S
    @eval function SymbolicUtils.promote_symtype(
        ::$(typeof(f)), T::Type{<:QField}, S::Type{<:QField}
    )
        return promote_type(T, S)
    end
end

SymbolicUtils.symtype(x::T) where {T<:QField} = T

Base.one(::T) where {T<:QField} = one(T)
Base.one(::Type{<:QField}) = 1
Base.isone(::QField) = false
Base.zero(::T) where {T<:QField} = zero(T)
Base.zero(::Type{<:QField}) = 0
Base.iszero(::QField) = false

#########################
# Field enum properties
#########################

name(ϕ::QSym) = ϕ.name

"""
    Regularisation `Plus` `Zero` `Minus`

Regularisation enum for the Keldysh quantum field. The regularisation is used to perform tadpole regularisation during the Wick contraction.
"""
@enum Regularisation begin
    Plus = 1
    Zero = 0
    Minus = -1
end
subtraction(x::NTuple{2,Regularisation}) = -(Int.(x)...)
function subtraction(x::Vector{Regularisation})
    @assert length(x) == 2
    return subtraction(Tuple(x))
end

"""
    KeldyshContour `Quantum` `Classical`

Keldysh contour enum for the Keldysh quantum field. The Keldysh contour is used to distinguish the field on quantum and classical contour.
"""
@enum KeldyshContour begin
    Quantum = 0
    Classical = 1
end
is_quantum(x::QSym) = iszero(Int(contour(x)))
is_classical(x::QSym) = isone(Int(contour(x)))

#########################
#       Position
#########################

"""
$(TYPEDEF)

Position type for the Keldysh quantum field.
The position is used to determine the coordinate of the field during the wick contraction.

Position has thee cases:
- [`In`](@ref): Position.index will be `typemax(Int8)`, which is used to represent the input field.Position
- [`Out`](@ref): Position.index will be `typemin(Int8)`, which is used to represent the output field.Position
- [`Bulk`](@ref): Position.index will be a positive integer
"""
struct Position
    index::Int8
end
"""
    Bulk(i::Int=1)
Create a `Bulk` position with a positive integer index `i`.
"""
function Bulk(i::Int=1)
    @assert i > 0 "Bulk index must be positive"
    return Position(i)
end
"""
    In()
Create a `Position` representing the input field.
"""
In() = Position(typemax(Int8))
"""
    Out()
Create a `Position` representing the output field.
"""
Out() = Position(typemin(Int8))

Base.isless(x::Position, y::Position) = x.index < y.index

is_bulk(x::Position) = typemin(Int8) < index(x) < typemax(Int8)
is_in(x::Position) = isequal(index(x), typemax(Int8))
is_out(x::Position) = isequal(index(x), typemin(Int8))

index(p::Position) = p.index
Base.Int(p::Position) = index(p)

# TODO make has_in and replace `In() ∈ positions`
