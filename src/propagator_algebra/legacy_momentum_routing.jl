# (* Legacy momentum routing compatibility *)

function _legacy_momentum_transform(gf::DressedPropagator{C,S,O,E1,E2}) where {C,S,O,E1,E2}
    keldysh = construct_momenta_from_gf(gf.keldysh)
    retarded = construct_momenta_from_gf(gf.retarded)
    advanced = construct_momenta_from_gf(gf.advanced)
    return DressedPropagator{C,S,O,E1,E2}(keldysh, retarded, advanced, gf.parameter)
end

function _legacy_momentum_transform(se::SelfEnergy{C,S,O,E1,E2}) where {C,S,O,E1,E2}
    keldysh = construct_momenta_from_self_energy(se.keldysh)
    retarded = construct_momenta_from_self_energy(se.retarded)
    advanced = construct_momenta_from_self_energy(se.advanced)
    return SelfEnergy{C,S,O,E1,E2}(keldysh, retarded, advanced, se.parameter)
end

function construct_momenta_from_gf(d::Diagram{S,E1,E2}) where {S,E1,E2}
    if iszero(E2)
        momenta = FixedVector{E1,Momenta}(
            Momenta(has_in(ps) || has_out(ps) ? 0 : 1) for ps in positions.(d.contractions)
        )
    else
        A = construct_linear_system(d.contractions)
        dep_idx, free_idx, P = solve_linear_system(A)
        momenta_dynamic = construct_momenta(dep_idx, free_idx, P)
        pushfirst!(momenta_dynamic, Momenta(0))
        push!(momenta_dynamic, Momenta(0))
        momenta = FixedVector{E1,Momenta}(momenta_dynamic)
    end
    return Diagram(d, momenta)
end
function construct_momenta_from_gf(d::Diagrams{C,S,E1,E2}) where {C,S,E1,E2}
    new_diagrams = Diagrams{C,S,E1,E2}()
    for (diagram, prefactor) in d
        new_diagram = construct_momenta_from_gf(diagram)
        push!(new_diagrams, new_diagram, prefactor)
    end
    return new_diagrams
end

function construct_momenta_from_self_energy(d::Diagram{S,E1,E2}) where {S,E1,E2}
    if iszero(E2)
        momenta = FixedVector{E1,Momenta}(
            Momenta(has_in(ps) || has_out(ps) ? 0 : 1) for ps in positions.(d.contractions)
        )
    else
        A = construct_linear_system(d.contractions)
        A = hcat([-1; 0], A)
        if !iseven(first(d.topology))
            A = hcat(A, [0; 1])
        else
            A = hcat(A, [1; 0])
        end

        dep_idx, free_idx, P = solve_linear_system(A)
        momenta_dynamic = construct_momenta(dep_idx, free_idx, P)
        momenta = FixedVector{E1,Momenta}(momenta_dynamic)
    end
    return Diagram(d, momenta)
end
function construct_momenta_from_self_energy(d::Diagrams{C,S,E1,E2}) where {C,S,E1,E2}
    new_diagrams = Diagrams{C,S,E1,E2}()
    for (diagram, prefactor) in d
        new_diagram = construct_momenta_from_self_energy(diagram)
        push!(new_diagrams, new_diagram, prefactor)
    end
    return new_diagrams
end

function construct_momenta(dep_idx, free_idx, P)
    l = size(P, 2)

    if isone(length(unique(eachrow(P))))
        return Momenta[Momenta(0) for _ in 1:l]
    end
    out = Vector{Momenta}(undef, l)

    idxs = free_idx .- 1
    idxs[end] = 0

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

function construct_linear_system(contractions)
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
    idxs = Int[1, idx]
    idxs_diff = setdiff(1:size(A, 2), idxs)
    A1 = A[:, idxs]
    A2 = A[:, idxs_diff]

    m, n = size(A1)
    if m != n
        error("A[:, dep_idx] must be square.")
    end
    if LinearAlgebra.det(A1) == 0
        error("A[:, dep_idx] is singular and not invertible.")
    end

    P = -(A1 \ A2)
    return idxs, idxs_diff, P
end
