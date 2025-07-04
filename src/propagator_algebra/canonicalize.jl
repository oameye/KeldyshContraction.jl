"""
    advanced_to_retarded(x::T) where {T<:SymbolicUtils.Symbolic}

Apply the transformation to change the advanced propagator to retarded:

``G^A(y, y)=-G^R(y, y)``

with ``y =(\\vec{y},t)``.
Note the expression is only valid for equal space-time coordinates.
"""
function advanced_to_retarded(
    x::Vector{Contraction}, prefactor::Number
)::Tuple{Vector{Contraction},Number}
    ff(x::Contraction) = is_advanced(x) && same_position(x)
    adv_idx = findall(ff, x)
    if isempty(adv_idx)
        return x, prefactor
    end
    x′ = deepcopy(x)
    # TODO or make x immutable or change in place
    for i in adv_idx
        prefactor *= -1
        x′[i] = adjoint(x[i])
    end
    return x′, prefactor
end

function sort_by_position_and_type(p::Contraction)::Float64
    if is_out(p)
        return -Inf
    elseif is_in(p)
        return Inf
    else
        i, j = integer_positions(p)
        type = Int(propagator_type(p...))
        return float(pairing(i, j) * 3 + type)
    end
end
sort_by_position_and_type(p::Edge)::Float64 = sort_by_position_and_type(fields(p))

function canonicalize(vs::Vector{Contraction})
    # TODO: assumes lower then third order
    canonical = is_canonical(vs)
    if !canonical
        vs = map(vs) do c
            map(c) do ψ
                if is_bulk(ψ)
                    idx = index(position(ψ))
                    ψ(Bulk(mod1(idx + 1, 2)))
                else
                    ψ
                end
            end
        end
    end
    return vs
end

"""
A canonical contraction (diagram) is one where the output field is connected to the
bulk position with the lowest index, canonically `index=1`.
TODO: How should higher orders be handled?
"""
function is_canonical(vs::Vector{Contraction})
    # TODO: assumes lower then third order
    idx_out = findfirst(c -> has_out(position.(c)), vs)
    if isnothing(idx_out)
        return true
    else
        return Bulk() ∈ position.(vs[idx_out]) ? true : false
    end
end
