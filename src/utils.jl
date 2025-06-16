sort_tuple(t) = t[1] < t[2] ? t : (t[2], t[1])

"""
    pairing(i,j)

The Hopcroft-Ullman pairing function, which bijectively maps two positive integers to a positive integer.

See https://mathworld.wolfram.com/PairingFunction.html
"""
@inline pairing(i::Integer, j::Integer) = div((i+j-2)*(i+j-1), 2) + i

make_real(x::Number) = SymbolicUtils._isreal(x) ? real(x) : x

_simplify(x::Complex) =
    if iszero(x.im)
        complex(x.re)
    elseif iszero(x.re)
        0.0+im*x.im
    else
        x
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
