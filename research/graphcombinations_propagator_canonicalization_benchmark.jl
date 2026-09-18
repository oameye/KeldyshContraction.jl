include(joinpath(@__DIR__, "graphcombinations_propagator_canonicalization_oracle.jl"))

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

function best_batch_time(f, repetitions; samples=7)
    best = Inf
    for _ in 1:samples
        elapsed = @elapsed for _ in 1:repetitions
            f()
        end
        best = min(best, elapsed)
    end
    return best / repetitions
end

function allocation_per_call(f, repetitions)
    bytes = @allocated for _ in 1:repetitions
        f()
    end
    return bytes / repetitions
end

@qfields ψ::Boson
ψc, ψq = ψ[KC.Classical], ψ[KC.Quantum]
@qfields η::Boson
ηc, ηq = η[KC.Classical], η[KC.Quantum]

benchmarks = Pair{String,Any}[]

push!(benchmarks, "2-node" => [
    KC.Contraction(ψc(KC.Out()), KC.bar(ψq)(KC.Bulk(1))),
    KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(2))),
    KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.In())),
])

push!(benchmarks, "3-node" => [
    KC.Contraction(ψc(KC.Out()), KC.bar(ψq)(KC.Bulk(2))),
    KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.Bulk(1))),
    KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.Bulk(3))),
    KC.Contraction(ψc(KC.Bulk(3)), KC.bar(ψq)(KC.Bulk(2))),
    KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(3))),
    KC.Contraction(ψc(KC.Bulk(3)), KC.bar(ψq)(KC.Bulk(1))),
    KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.In())),
])

push!(benchmarks, "4-node" => [
    KC.Contraction(ψc(KC.Out()), KC.bar(ψq)(KC.Bulk(3))),
    KC.Contraction(ψc(KC.Bulk(3)), KC.bar(ψq)(KC.Bulk(1))),
    KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(4))),
    KC.Contraction(ψc(KC.Bulk(4)), KC.bar(ψq)(KC.Bulk(2))),
    KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.Bulk(3))),
    KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.In())),
])

push!(benchmarks, "5-node" => [
    KC.Contraction(ψc(KC.Out()), KC.bar(ψq)(KC.Bulk(5))),
    KC.Contraction(ψc(KC.Bulk(5)), KC.bar(ψq)(KC.Bulk(1))),
    KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(2))),
    KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.Bulk(3))),
    KC.Contraction(ψc(KC.Bulk(3)), KC.bar(ψq)(KC.Bulk(4))),
    KC.Contraction(ψc(KC.Bulk(4)), KC.bar(ψq)(KC.Bulk(5))),
    KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(3))),
    KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.Bulk(4))),
    KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.In())),
])

push!(benchmarks, "colored-2-node" => [
    KC.Contraction(ψc(KC.Out()), KC.bar(ψq)(KC.Bulk(1))),
    KC.Contraction(ηc(KC.Bulk(1)), KC.bar(ηq)(KC.Bulk(2))),
    KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.In())),
])

push!(benchmarks, "parallel-2-node" => [
    KC.Contraction(ψc(KC.Out()), KC.bar(ψq)(KC.Bulk(1))),
    KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(2))),
    KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(2))),
    KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.In())),
])

capacity = maximum(length(KC.canonicalization_positions(vs)) + length(vs) for (_, vs) in benchmarks)
scratch = GCPhysicalCanonicalizationWorkspace(capacity)

println()
println("production-shaped warmed canonicalization benchmark")
println("  reusable GC capacity: ", capacity)
println("  direct graph storage and witness transport are caller-reused")
println("  times are best-of-7 batch means; allocations are warmed batch means")

for (label, vs) in benchmarks
    gc_call = () -> gc_physical_canonicalize!(scratch, vs)
    nauty_call = () -> KC.canonicalize(vs)

    for _ in 1:100
        gc_call()
        nauty_call()
    end

    repetitions = label == "5-node" ? 200 : 1000
    gc_time = best_batch_time(gc_call, repetitions)
    nauty_time = best_batch_time(nauty_call, repetitions)
    gc_alloc = allocation_per_call(gc_call, repetitions)
    nauty_alloc = allocation_per_call(nauty_call, repetitions)

    println(
        "  ", label,
        ": GC/Nauty time=", round(gc_time / nauty_time; digits=3), "x",
        " alloc=", round(gc_alloc / nauty_alloc; digits=3), "x",
        " [GC ", round(gc_time * 1e6; digits=2), " us, ", round(gc_alloc; digits=1), " B",
        "; Nauty ", round(nauty_time * 1e6; digits=2), " us, ", round(nauty_alloc; digits=1), " B]",
    )
end
