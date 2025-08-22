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

function canonicalize(vs::Vector{Contraction}) # TODO: cache this on top of Vector{Ntuple{2,Position}}
    all_positions = vcat([collect(position.(c)) for c in vs]...)
    bulk_nodes = unique(filter(is_bulk, all_positions))
    isempty(bulk_nodes) && return vs

    # Find which bulk node is connected to Out()
    out_contraction = findfirst(c -> any(pos -> pos == Out(), position.(c)), vs)
    if isnothing(out_contraction)
        return vs  # No Out() connection found
    end

    out_positions = position.(vs[out_contraction])
    bulk_positions = filter(is_bulk, out_positions)
    if isempty(bulk_positions)
        return vs  # Out() not connected to bulk
    end

    out_connected_bulk = first(bulk_positions)

    # Only consider permutations where the bulk node connected to Out() becomes Bulk(1)
    best_perm = nothing
    best_sig = nothing

    for perm in SmallCollections.Combinatorics.permutations(length(bulk_nodes))
        # Check if this permutation puts the Out()-connected bulk node at position 1
        out_bulk_idx = findfirst(==(out_connected_bulk), bulk_nodes)
        if perm[out_bulk_idx] == 1
            sig = graph_signature(vs, bulk_nodes, perm)
            if best_sig === nothing || sig < best_sig
                best_sig = sig
                best_perm = perm
            end
        end
    end

    return best_perm === nothing ? vs : apply_permutation(vs, bulk_nodes, best_perm)
end

function graph_signature(vs::Vector{Contraction}, bulk_nodes, perm)
    # Build adjacency matrix with ordered nodes: Out, Bulk(perm[1]), Bulk(perm[2]), ..., In
    node_map = Dict{Position,Int}()
    node_map[Out()] = 1
    for (i, bulk_idx) in enumerate(perm)
        node_map[bulk_nodes[bulk_idx]] = i + 1
    end
    node_map[In()] = length(bulk_nodes) + 2

    n = length(bulk_nodes) + 2
    adj = zeros(Int, n, n)

    for c in vs
        pos1, pos2 = position.(c)
        i, j = node_map[pos1], node_map[pos2]
        adj[i, j] += 1
        adj[j, i] += 1  # undirected
    end

    return vec(adj)  # flatten to vector for comparison
end

function apply_permutation(vs::Vector{Contraction}, bulk_nodes, perm)
    # Create mapping from old bulk nodes to new canonical numbering
    perm_map = Dict{Position,Position}()
    for (i, bulk_idx) in enumerate(perm)
        perm_map[bulk_nodes[bulk_idx]] = Bulk(i)
    end

    return map(vs) do c
        map(c) do ψ
            pos = position(ψ)
            if is_bulk(pos) && haskey(perm_map, pos)
                ψ(perm_map[pos])
            else
                ψ
            end
        end
    end
end

function is_canonical(vs::Vector{Contraction})
    isequal(canonicalize(vs), vs)
end
