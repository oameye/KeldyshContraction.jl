sort_tuple(t) = t[1] < t[2] ? t : (t[2], t[1])

"""
    pairing(i,j)

The Hopcroft-Ullman pairing function, which bijectively maps two positive integers to a positive integer.

See https://mathworld.wolfram.com/PairingFunction.html
"""
@inline pairing(i::Int, j::Int)::Int = div((i+j-2)*(i+j-1), 2) + i

make_real(x::Number) = SymbolicUtils._isreal(x) ? real(x) : x

_simplify(x::Complex) =
    if iszero(x.im)
        complex(x.re)
    elseif iszero(x.re)
        0.0+im*x.im
    else
        x
    end

function expand_coefficients(v::Vector{Int})
    v′ = deepcopy(v)
    o = Vector{Int}()
    while sum(v′) > 0
        for idx in eachindex(v′)
            if v′[idx] > 0
                push!(o, idx)
                v′[idx] -= 1
                break
            end
        end
    end
    return o
end # TODO: need better name
