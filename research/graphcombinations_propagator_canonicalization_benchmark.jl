include(joinpath(@__DIR__, "graphcombinations_propagator_canonicalization_oracle.jl"))

struct GCPhysicalCanonicalizationWorkspace
    search::GC.DirectedCanonicalizationWorkspace
    result::GC.DirectedCanonicalizationBuffer
end

function GCPhysicalCanonicalizationWorkspace(capacity::Integer)
    return GCPhysicalCanonicalizationWorkspace(
        GC.DirectedCanonicalizationWorkspace(capacity),
        GC.DirectedCanonicalizationBuffer(capacity),
    )
end

function gc_physical_permutation!(scratch::GCPhysicalCanonicalizationWorkspace, vs, graph_positions)
    direct_graph, direct_colors, simple = _gc_direct_graph(vs, graph_positions)
    GC.canonicalize_directed!(scratch.result, scratch.search, direct_graph, direct_colors)

    use_direct = isone(GC.canonical_automorphism_order(scratch.result)) ||
                 (simple && KC.uniform_coloring(vs))
    if !use_direct
        colored_graph, colored_colors = _gc_colored_graph(vs, graph_positions)
        GC.canonicalize_directed!(scratch.result, scratch.search, colored_graph, colored_colors)
    end

    return [GC.original_vertex(scratch.result, rank) for rank in 1:scratch.result.num_vertices]
end

function gc_physical_canonicalize!(scratch::GCPhysicalCanonicalizationWorkspace, vs::Vector{T}) where {T}
    isempty(vs) && return copy(vs)
    graph_positions = KC.canonicalization_positions(vs)
    permutation = gc_physical_permutation!(scratch, vs, graph_positions)
    mapping = KC.make_permutation_dict(permutation, graph_positions, vs)
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
