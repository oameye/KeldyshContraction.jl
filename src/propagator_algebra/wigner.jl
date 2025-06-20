# (* Experimental *)

function wigner_transform(gf::DressedPropagator)
    # Apply Wigner transform to each component
    keldysh = construct_momenta_from_gf(gf.keldysh)
    retarded = construct_momenta_from_gf(gf.retarded)
    advanced = construct_momenta_from_gf(gf.advanced)

    # Return a new DressedPropagator with the transformed components
    return DressedPropagator(keldysh, retarded, advanced, gf.order, gf.parameter)
end

function wigner_transform(se::SelfEnergy)
    # Apply Wigner transform to each component
    keldysh = construct_momenta_from_self_energy(se.keldysh)
    retarded = construct_momenta_from_self_energy(se.retarded)
    advanced = construct_momenta_from_self_energy(se.advanced)

    # Return a new DressedPropagator with the transformed components
    return SelfEnergy(keldysh, retarded, advanced, se.order, se.parameter)
end

function construct_momenta_from_gf(d::Diagram{E1,E2}) where {E1,E2}
    if iszero(E2)
        momenta = map(positions.(d.contractions)) do ps
            if has_in(ps) || has_out(ps)
                m = Momenta(0)
            else
                m = Momenta(1)
            end
            m
        end
    else
        A = construct_linear_system(d.contractions)
        dep_idx, free_idx, P = solve_linear_system(A)
        momenta = construct_momenta(dep_idx, free_idx, P)
        pushfirst!(momenta, Momenta(0))
        push!(momenta, Momenta(0))
        momenta = FixedVector(momenta)
    end
    return Diagram(d, momenta)
end
function construct_momenta_from_gf(d::Diagrams{E1,E2}) where {E1,E2}
    new_diagrams = Diagrams{E1,E2}()
    for (diagram, prefactor) in d
        new_diagram = construct_momenta_from_gf(diagram)
        push!(new_diagrams, new_diagram, prefactor)
    end
    return new_diagrams
end # TODO: can this be faster?

function construct_momenta_from_self_energy(d::Diagram{E1,E2}) where {E1,E2}
    if iszero(E2)
        momenta = map(positions.(d.contractions)) do ps
            if has_in(ps) || has_out(ps)
                m = Momenta(0)
            else
                m = Momenta(1)
            end
            m
        end
    else
        A = construct_linear_system(d.contractions)
        A = hcat([-1; 0], A) # canonicalize diagram
        if !iseven(first(d.topology))
            A = hcat(A, [0; 1]) #
        else
            A = hcat(A, [1; 0]) # canonicalize diagram
        end # TODO instead replace with topology -> positions function

        dep_idx, free_idx, P = solve_linear_system(A)
        momenta = construct_momenta(dep_idx, free_idx, P)
        momenta = FixedVector(momenta)
    end
    return Diagram(d, momenta)
end
function construct_momenta_from_self_energy(d::Diagrams{E1,E2}) where {E1,E2}
    new_diagrams = Diagrams{E1,E2}()
    for (diagram, prefactor) in d
        new_diagram = construct_momenta_from_self_energy(diagram)
        push!(new_diagrams, new_diagram, prefactor)
    end
    return new_diagrams
end # TODO: can this be faster?

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

# TODO topology -> positions
# function positions(topo::FixedVector{E,Int}) where {E}
#     return map(e -> e.position, topo)
# end
