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

function make_permutation_dict(perm, max_label, has_in)
    if has_in
        l = length(perm)
        tracker = 0
        last_index = findfirst(==(max_label), perm)
        first_index = findfirst(==(1), perm)
        dict = Dict{Position,Position}()
        for i in 1:l
            if i == first_index || i == last_index
                tracker += 1
            else
                dict[Bulk(perm[i] - 1)] = Bulk(i - tracker)
            end
        end
    else # vacuum diagram
        dict = Dict{Position,Position}(Bulk(perm[i]) => Bulk(i) for i in 1:length(perm))
    end

    return dict
end
function canonicalize(vs)
    graph, max_label, has_in = make_graph(NautyGraphs.NautyDiGraph, vs)
    perm = NautyGraphs.canonical_permutation(graph)
    permutation_map = make_permutation_dict(perm, max_label, has_in)

    canonical_vs = map(vs) do c
        map(c) do ψ
            pos = position(ψ)
            if is_bulk(pos) && haskey(permutation_map, pos)
                ψ(permutation_map[pos])
            else
                ψ
            end
        end
    end
    return canonical_vs
end
