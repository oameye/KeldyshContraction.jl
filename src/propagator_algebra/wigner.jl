struct Momentum
    index::Int8
    function Momentum(index)
        index < 0 && error("Momentum index must be non-negative, got $index")
        new(index)
    end
end

struct Momenta
    prefactors::Vector{Int}
    Momenta::Vector{Momentum}
    function Momenta(prefactors::Vector{Int}, Momenta::Vector{Momentum})
        @assert length(prefactors) == length(Momenta) "Length of prefactors and momenta must match"
        idxs = findall(x->!iszero(x), prefactors)
        new(prefactors[idxs], Momenta[idxs])
    end
end
function Momenta(idx::Int)
    Momenta([1], [Momentum(idx)])
end

function construct_momenta(dep_idx, free_idx, P)
    l = size(P, 2)

    if isone(length(unique(eachrow(P))))
        return [Momenta(0) for _ in 1:l]
    end
    out = Vector{Momenta}(undef, l)

    idxs = free_idx .- 1
    idxs[end] = 0 # last index is external momentum

    for idx in idxs
        if idx == 0
            continue
        end
        out[idx] = Momenta(idx)
    end
    ms = map(Momentum, idxs)
    out[last(dep_idx) - 1] = Momenta(Vector{Int}(P[2, :]), ms)
    return out
end

function construct_linear_system(contractions)::Matrix{Int}
    # TODO supports only second-order contractions
    A = zeros(Int, 2, length(contractions))
    pos = positions.(contractions)
    for (j, ps) in enumerate(pos)
        for (i, p) in enumerate(ps)
            idx = index(p)
            if idx == typemax(Int8) || idx == typemin(Int8)
                continue
            elseif !iszero(A[idx, j])
                A[idx, j] = 0
            else
                A[idx, j] = bool_to_index(isone(i))
            end
        end
    end
    return A
end

function solve_linear_system(A::Matrix{Int})
    # Extract submatrices A1 and A2
    l = size(A, 2)
    idx = l
    for diff_ in 1:size(A, 2)
        idx = l - diff_
        if iszero(A[:, idx])
            continue
        else
            break
        end
    end
    if idx == 1
        error("No two dependent index found in the linear system.")
    end
    idxs = [1, idx]
    idxs_diff = setdiff(1:size(A, 2), idxs)
    A1 = A[:, idxs]
    A2 = A[:, idxs_diff]

    # Check that A1 is square and invertible
    m, n = size(A1)
    if m != n
        error("A[:, dep_idx] must be square.")
    end
    if LinearAlgebra.det(A1) == 0
        error("A[:, dep_idx] is singular and not invertible.")
    end

    # Compute P = -inv(A1) * A2
    P = -inv(A1) * A2
    return idxs, idxs_diff, P
end

bool_to_index(x::Bool) = 2 * x - 1
