# KC-owned causal pruning for the GraphCombinations colored-port traversal.
#
# A completed causal cycle is extension-monotone, so it is safe to reject the raw child before
# GC canonicalizes it. Pair-local Keldysh semantics are compiled once into an opaque port-cell
# direction table; the hot policy works only on compact GC edges.

const _GC_NONCAUSAL = Int8(0)
const _GC_ADVANCED = Int8(1)
const _GC_RETARDED = Int8(2)

struct _GCCausalPortPolicy
    kinds::Array{Int8,4}
end

@inline function _gc_causal_kind(contraction::Contraction)::Int8
    propagator = propagator_type(contraction...)
    is_advanced(propagator) && return _GC_ADVANCED
    is_retarded(propagator) && return _GC_RETARDED
    return _GC_NONCAUSAL
end

function _gc_causal_policy(lookup::Dict{NTuple{4,Int},Contraction{S}}) where {S<:Statistics}
    isempty(lookup) && return _GCCausalPortPolicy(zeros(Int8, 0, 0, 0, 0))

    max_source_vertex = maximum(cell[1] for cell in keys(lookup))
    max_source_color = maximum(cell[2] for cell in keys(lookup))
    max_target_vertex = maximum(cell[3] for cell in keys(lookup))
    max_target_color = maximum(cell[4] for cell in keys(lookup))
    kinds = zeros(
        Int8,
        max_source_vertex,
        max_source_color,
        max_target_vertex,
        max_target_color,
    )
    for (cell, contraction) in lookup
        kinds[cell...] = _gc_causal_kind(contraction)
    end
    return _GCCausalPortPolicy(kinds)
end

@inline function _gc_edge_kind(policy::_GCCausalPortPolicy, edge::GC.ColoredPortEdge)::Int8
    return policy.kinds[edge.source, edge.source_color, edge.target, edge.target_color]
end

@inline function _gc_causal_vertices(edge::GC.ColoredPortEdge, kind::Int8)
    if kind == _GC_ADVANCED
        return edge.source, edge.target
    elseif kind == _GC_RETARDED
        return edge.target, edge.source
    end
    return 0, 0
end

@inline _gc_port_vertex_bit(vertex::Int) = one(UInt128) << (vertex - 1)

function _gc_has_causal_path(
    policy::_GCCausalPortPolicy, edges, source::Int, target::Int, edge_count::Int
)::Bool
    reachable = _gc_port_vertex_bit(source)
    target_bit = _gc_port_vertex_bit(target)
    !iszero(reachable & target_bit) && return true

    previous = zero(UInt128)
    while reachable != previous
        previous = reachable
        @inbounds for i in 1:edge_count
            edge = edges[i]
            kind = _gc_edge_kind(policy, edge)
            kind == _GC_NONCAUSAL && continue
            causal_source, causal_target = _gc_causal_vertices(edge, kind)
            causal_source == causal_target && continue
            if !iszero(reachable & _gc_port_vertex_bit(causal_source))
                reachable |= _gc_port_vertex_bit(causal_target)
            end
        end
        !iszero(reachable & target_bit) && return true
    end
    return false
end

function (policy::_GCCausalPortPolicy)(state::GC.ColoredPortState)::Bool
    edges = GC.port_edges(state)
    isempty(edges) && return true

    added = last(edges)
    kind = _gc_edge_kind(policy, added)
    kind == _GC_NONCAUSAL && return true

    source, target = _gc_causal_vertices(added, kind)
    parent_edge_count = length(edges) - 1
    if source == target
        @inbounds for i in 1:parent_edge_count
            edge = edges[i]
            previous_kind = _gc_edge_kind(policy, edge)
            previous_kind == _GC_NONCAUSAL && continue
            previous_source, previous_target = _gc_causal_vertices(edge, previous_kind)
            previous_source == source == previous_target && return false
        end
        return true
    end

    return !_gc_has_causal_path(policy, edges, target, source, parent_edge_count)
end

function _gc_wick_contraction_with_stats(
    args_nc::Vector{Field{S}},
    ::Val{E},
    ::Val{E2};
    regularise=true,
    _set_reg_to_zero=false,
    simplify=false,
    causal_pruning=true,
) where {S<:Union{Boson,Fermion},E,E2}
    problem, lookup, transport = _gc_build_wick_problem(
        args_nc, Val(E); regularise, _set_reg_to_zero
    )
    policy = causal_pruning ? _gc_causal_policy(lookup) : GC.AcceptAllPortPolicy()
    completions, stats = GC.generate_weighted_with_stats(problem; policy, transport)

    canonical_scratch = physical_canonicalization_workspace(Val(E))
    canonical_weights = Dict{FixedVector{E,Contraction{S}},BigInt}()
    for completion in completions
        contractions = _gc_completion_contractions(completion, lookup, Val(E))
        if causal_pruning
            is_connected(contractions) || continue
        else
            passes_wick_filters(contractions) || continue
        end

        final_contractions, simplification_sign = if simplify
            advanced_to_retarded(contractions, 1)
        else
            contractions, 1
        end
        canonical = canonicalize(final_contractions, canonical_scratch)
        canonical_key = sorted_wick_key(canonical, Val(E))
        final_weight = completion.weight * Int(simplification_sign)
        canonical_weights[canonical_key] =
            get(canonical_weights, canonical_key, big(0)) + final_weight
    end

    wick_pairings = Tuple{WickPairing{S,E},FixedVector{E2,Int},Int}[]
    sizehint!(wick_pairings, length(canonical_weights))
    for (canonical, weight) in canonical_weights
        iszero(weight) && continue
        contractions = Contraction{S}[contraction for contraction in canonical]
        topology = legacy_topology(contractions, Val(E2))
        pairing = WickPairing(contractions, Int8(sign(weight)), Val(E))
        push!(wick_pairings, (pairing, topology, Int(abs(weight))))
    end
    return wick_pairings, stats
end
