sort_tuple(t) = t[1] < t[2] ? t : (t[2], t[1])

"""
    pairing(i,j)

The Hopcroft-Ullman pairing function, which bijectively maps two positive integers to a positive integer.

See https://mathworld.wolfram.com/PairingFunction.html
"""
@inline pairing(i::Integer, j::Integer) = div((i+j-2)*(i+j-1), 2) + i

"""
maps edge (i,j) to a unique integer in range 1:max_edges(n).
It is different from Hopcroft-Ullman pairing as it assumes that i ≠ j and the number of
edges is used.
"""
@inline function edge_to_index(i::Integer, j::Integer, n::Integer)::Integer
    @boundscheck begin
        1 ≤ i ≤ n || throw(ArgumentError("i must be in 1:$n, got $i"))
        1 ≤ j ≤ n || throw(ArgumentError("j must be in 1:$n, got $j"))
        i ≠ j || throw(ArgumentError("i and j must be different"))
    end
    i, j = minmax(i, j)                 # ensure i < j
    return (i - 1) * n - ((i - 1) * i) ÷ 2 + (j - i)
end

bool_to_index(x::Bool) = 2 * x - 1

# Numeric coefficients are part of the canonical symbolic representation.  In
# particular, IEEE signed zero must not leak into `isequal`/`hash` semantics:
# `-0.0` and `0.0` are algebraically identical even though `isequal` distinguishes
# their bit patterns.  Normalize zero components while preserving the concrete
# coefficient type.
_simplify(x::Number) = x
_simplify(x::T) where {T<:AbstractFloat} = iszero(x) ? zero(T) : x
function _simplify(x::Complex{T}) where {T}
    re = iszero(x.re) ? zero(x.re) : x.re
    im = iszero(x.im) ? zero(x.im) : x.im
    return complex(re, im)
end

"""
    indices_from_counts(counts::Vector{Int})

Convert a vector of counts into a flat vector of indices where each index appears
according to its count.

## Examples
```jldoctest
julia> using KeldyshContraction: indices_from_counts

julia> indices_from_counts([2, 0, 1])
3-element Vector{Int64}:
 1
 1
 3

julia> indices_from_counts([1, 2, 1])
4-element Vector{Int64}:
 1
 2
 2
 3
```

## Notes
This function is commonly used to expand multinomial coefficients from
`Combinatorics.multiexponents` into index lists for perturbation theory calculations.
"""
function indices_from_counts(counts::Vector{Int})
    counts′ = deepcopy(counts)
    indices = Vector{Int}()
    while sum(counts′) > 0
        for idx in eachindex(counts′)
            if counts′[idx] > 0
                push!(indices, idx)
                counts′[idx] -= 1
                break
            end
        end
    end
    return indices
end
