import GraphCombinations as GC
import KeldyshContraction as KC

struct GCNativePhysicalCanonicalizationWorkspace
    direct::GCPhysicalCanonicalizationWorkspace
    relation_search::GC.DirectedRelationCanonicalizationWorkspace
    relation_result::GC.DirectedRelationCanonicalizationBuffer
    relation_graph::GC.DirectedRelationGraphBuffer
end

function GCNativePhysicalCanonicalizationWorkspace(capacity::Integer)
    n = Int(capacity)
    return GCNativePhysicalCanonicalizationWorkspace(
        GCPhysicalCanonicalizationWorkspace(n),
        GC.DirectedRelationCanonicalizationWorkspace(n, n),
        GC.DirectedRelationCanonicalizationBuffer(n, n; materialize_canonical=false),
        GC.DirectedRelationGraphBuffer(n, n),
    )
end

function gc_native_relation_graph!(
    scratch::GCNativePhysicalCanonicalizationWorkspace, vs, graph_positions
)
    colors = KC.propagator_colors(vs)
    n = length(graph_positions)
    nr = length(colors)
    graph = scratch.relation_graph
    n <= graph.vertex_capacity ||
        throw(DimensionMismatch("native relation vertex capacity is too small"))
    nr <= graph.relation_capacity ||
        throw(DimensionMismatch("native relation color capacity is too small"))

    active = nr * n * n
    @inbounds for slot in 1:active
        graph.multiplicities[slot] = 0
    end
    graph.num_vertices = n
    graph.num_relations = nr

    @inbounds for item in vs
        relation = searchsortedfirst(colors, KC.propagator_color(item))
        source_position, target_position = KC.positions(item)
        source = KC.position_vertex(graph_positions, source_position)
        target = KC.position_vertex(graph_positions, target_position)
        slot = ((relation - 1) * n + source - 1) * n + target
        graph.multiplicities[slot] += 1
    end
    return graph
end

function gc_native_relation_witness!(
    scratch::GCNativePhysicalCanonicalizationWorkspace,
    vs,
    graph_positions,
    vertex_colors=KC.position_labels(graph_positions),
)
    graph = gc_native_relation_graph!(scratch, vs, graph_positions)
    GC.canonicalize_directed_relations!(
        scratch.relation_result, scratch.relation_search, graph, vertex_colors
    )
    return scratch.relation_result
end

function gc_native_physical_canonicalize!(
    scratch::GCNativePhysicalCanonicalizationWorkspace, vs::Vector{T}
) where {T}
    isempty(vs) && return copy(vs)
    graph_positions = KC.canonicalization_positions(vs)
    witness = gc_native_relation_witness!(scratch, vs, graph_positions)
    mapping = gc_make_permutation_dict(witness, graph_positions, vs)
    return T[KC.relabel_bulk_positions(item, mapping) for item in vs]
end

function gc_native_hybrid_physical_canonicalize!(
    scratch::GCNativePhysicalCanonicalizationWorkspace, vs::Vector{T}
) where {T}
    isempty(vs) && return copy(vs)
    graph_positions = KC.canonicalization_positions(vs)
    direct_graph, vertex_colors, simple = gc_direct_graph!(
        scratch.direct, vs, graph_positions
    )
    GC.canonicalize_directed!(
        scratch.direct.result, scratch.direct.search, direct_graph, vertex_colors
    )

    use_direct =
        isone(GC.canonical_automorphism_order(scratch.direct.result)) ||
        (simple && KC.uniform_coloring(vs))
    if use_direct
        mapping = gc_make_permutation_dict(scratch.direct.result, graph_positions, vs)
        return T[KC.relabel_bulk_positions(item, mapping) for item in vs]
    end

    witness = gc_native_relation_witness!(scratch, vs, graph_positions, vertex_colors)
    mapping = gc_make_permutation_dict(witness, graph_positions, vs)
    return T[KC.relabel_bulk_positions(item, mapping) for item in vs]
end
