import GraphCombinations as GC
import KeldyshContraction as KC

struct GCPhysicalCanonicalizationWorkspace
    search::GC.DirectedCanonicalizationWorkspace
    result::GC.DirectedCanonicalizationBuffer
    multiplicities::Vector{Int}
end

function GCPhysicalCanonicalizationWorkspace(capacity::Integer)
    n = Int(capacity)
    return GCPhysicalCanonicalizationWorkspace(
        GC.DirectedCanonicalizationWorkspace(n),
        GC.DirectedCanonicalizationBuffer(n),
        zeros(Int, n * n),
    )
end

@inline function clear_active_multiplicities!(scratch::GCPhysicalCanonicalizationWorkspace, n::Int)
    @inbounds for i in 1:(n * n)
        scratch.multiplicities[i] = 0
    end
    return nothing
end

function gc_direct_graph!(scratch::GCPhysicalCanonicalizationWorkspace, vs, graph_positions)
    n = length(graph_positions)
    clear_active_multiplicities!(scratch, n)
    simple = true
    @inbounds for item in vs
        source_position, target_position = KC.positions(item)
        source = KC.position_vertex(graph_positions, source_position)
        target = KC.position_vertex(graph_positions, target_position)
        slot = (source - 1) * n + target
        simple &= iszero(scratch.multiplicities[slot])
        scratch.multiplicities[slot] += 1
    end
    graph = GC.DirectedGCGraph(n, scratch.multiplicities)
    return graph, KC.position_labels(graph_positions), simple
end

function gc_colored_graph!(scratch::GCPhysicalCanonicalizationWorkspace, vs, graph_positions)
    colors = KC.propagator_colors(vs)
    npositions = length(graph_positions)
    n = npositions + length(vs)
    clear_active_multiplicities!(scratch, n)

    labels = Vector{Int}(undef, n)
    position_colors = KC.position_labels(graph_positions)
    copyto!(labels, 1, position_colors, 1, npositions)

    @inbounds for (i, item) in enumerate(vs)
        color_index = searchsortedfirst(colors, KC.propagator_color(item))
        edge_vertex = npositions + i
        labels[edge_vertex] = 3 + color_index
        source_position, target_position = KC.positions(item)
        source = KC.position_vertex(graph_positions, source_position)
        target = KC.position_vertex(graph_positions, target_position)
        scratch.multiplicities[(source - 1) * n + edge_vertex] += 1
        scratch.multiplicities[(edge_vertex - 1) * n + target] += 1
    end
    return GC.DirectedGCGraph(n, scratch.multiplicities), labels
end

function gc_physical_witness!(scratch::GCPhysicalCanonicalizationWorkspace, vs, graph_positions)
    direct_graph, direct_colors, simple = gc_direct_graph!(scratch, vs, graph_positions)
    GC.canonicalize_directed!(scratch.result, scratch.search, direct_graph, direct_colors)

    use_direct = isone(GC.canonical_automorphism_order(scratch.result)) ||
                 (simple && KC.uniform_coloring(vs))
    if !use_direct
        colored_graph, colored_colors = gc_colored_graph!(scratch, vs, graph_positions)
        GC.canonicalize_directed!(scratch.result, scratch.search, colored_graph, colored_colors)
    end
    return scratch.result
end

function gc_make_permutation_dict(buffer, graph_positions, vs)
    npositions = length(graph_positions)
    canonical_bulk = KC.Position[]
    sizehint!(canonical_bulk, npositions)
    @inbounds for rank in 1:buffer.num_vertices
        original_vertex = GC.original_vertex(buffer, rank)
        original_vertex <= npositions || continue
        old_position = graph_positions[original_vertex]
        KC.is_bulk(old_position) && push!(canonical_bulk, old_position)
    end

    anchors = KC.out_bulk_positions(vs)
    mapping = Dict{KC.Position,KC.Position}()
    sizehint!(mapping, length(canonical_bulk))
    bulk_index = 0
    for old_position in canonical_bulk
        old_position in anchors || continue
        bulk_index += 1
        mapping[old_position] = KC.Bulk(bulk_index)
    end
    for old_position in canonical_bulk
        old_position in anchors && continue
        bulk_index += 1
        mapping[old_position] = KC.Bulk(bulk_index)
    end
    return mapping
end

function gc_physical_canonicalize!(scratch::GCPhysicalCanonicalizationWorkspace, vs::Vector{T}) where {T}
    isempty(vs) && return copy(vs)
    graph_positions = KC.canonicalization_positions(vs)
    witness = gc_physical_witness!(scratch, vs, graph_positions)
    mapping = gc_make_permutation_dict(witness, graph_positions, vs)
    return T[KC.relabel_bulk_positions(item, mapping) for item in vs]
end
