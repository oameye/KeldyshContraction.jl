const DERIVATIVE_AXES = (:t, :x, :y, :z)
const DERIVATIVE_ORDER_BITS = 8
const DERIVATIVE_ORDER_MASK = UInt32(0xff)

"""Canonical coordinate-derivative multi-index carried by a `Field`."""
struct DerivativeMultiIndex
    orders::UInt32
end

DerivativeMultiIndex() = DerivativeMultiIndex(UInt32(0))

@inline function derivative_axis_slot(axis::Symbol)::Int
    axis === :t && return 1
    axis === :x && return 2
    axis === :y && return 3
    axis === :z && return 4
    throw(
        ArgumentError("unsupported derivative axis: $axis; expected one of :t, :x, :y, :z")
    )
end

@inline derivative_shift(slot::Int) = DERIVATIVE_ORDER_BITS * (slot - 1)

@inline function derivative_order(d::DerivativeMultiIndex, slot::Int)::UInt8
    shift = derivative_shift(slot)
    return UInt8((d.orders >> shift) & DERIVATIVE_ORDER_MASK)
end

function DerivativeMultiIndex(axes::AbstractVector{Symbol})
    result = DerivativeMultiIndex()
    for axis in axes
        result = append_derivative(result, axis)
    end
    return result
end

function Base.length(d::DerivativeMultiIndex)
    total = 0
    @inbounds for slot in eachindex(DERIVATIVE_AXES)
        total += Int(derivative_order(d, slot))
    end
    return total
end
@inline Base.isempty(d::DerivativeMultiIndex) = iszero(d.orders)

function Base.getindex(d::DerivativeMultiIndex, i::Int)
    1 <= i <= length(d) || throw(BoundsError(d, i))
    seen = 0
    @inbounds for slot in eachindex(DERIVATIVE_AXES)
        seen += Int(derivative_order(d, slot))
        i <= seen && return DERIVATIVE_AXES[slot]
    end
    return throw(BoundsError(d, i))
end

function Base.iterate(d::DerivativeMultiIndex, state::Int=1)
    state > length(d) && return nothing
    return d[state], state + 1
end
Base.eltype(::Type{DerivativeMultiIndex}) = Symbol
Base.IteratorSize(::Type{DerivativeMultiIndex}) = Base.HasLength()

@inline Base.isequal(a::DerivativeMultiIndex, b::DerivativeMultiIndex) =
    a.orders == b.orders
Base.:(==)(a::DerivativeMultiIndex, b::DerivativeMultiIndex) = isequal(a, b)
@inline Base.hash(d::DerivativeMultiIndex, h::UInt) =
    hash(DerivativeMultiIndex, hash(d.orders, h))

@inline function Base.isless(a::DerivativeMultiIndex, b::DerivativeMultiIndex)
    a.orders == b.orders && return false
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        ai = a[i]
        bi = b[i]
        ai === bi || return isless(ai, bi)
    end
    return length(a) < length(b)
end

@inline function append_derivative(d::DerivativeMultiIndex, axis::Symbol)
    slot = derivative_axis_slot(axis)
    order = derivative_order(d, slot)
    order == typemax(UInt8) &&
        throw(OverflowError("derivative order exceeds $(typemax(UInt8)) along $axis"))
    increment = UInt32(1) << derivative_shift(slot)
    return DerivativeMultiIndex(d.orders + increment)
end

"""Return an owned canonical derivative-axis sequence."""
function derivative_axes(d::DerivativeMultiIndex)
    axes = Symbol[]
    sizehint!(axes, length(d))
    @inbounds for slot in eachindex(DERIVATIVE_AXES)
        for _ in 1:Int(derivative_order(d, slot))
            push!(axes, DERIVATIVE_AXES[slot])
        end
    end
    return axes
end
