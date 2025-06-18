
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
