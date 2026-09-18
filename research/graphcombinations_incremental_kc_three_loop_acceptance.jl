using BenchmarkTools
using Test

import GraphCombinations as GC
import KeldyshContraction as KC

include(ENV["KC_CORPUS_DEFS"])
include(ENV["KC_MATRIXFREE_DEFS"])

const ACCEPTANCE_LOOP_COUNT = 3

function _incremental_stats(workspace)
    candidate = workspace.canonicalization.workspace
    base = candidate.orbit.base
    return (
        vertices=workspace.canonicalization.graph.num_vertices,
        levels=base.levels,
        generated=base.generated_nodes,
        retained=base.retained_nodes,
        discarded_by_trace=base.discarded_by_trace,
        max_frontier=base.maximum_frontier,
        experimental_paths=base.experimental_paths,
        image_matches=base.exact_image_matches,
        quotient_discards=base.quotient_discards,
        generators=candidate.orbit.generators,
        orbit_skips=candidate.orbit.orbit_path_skips,
    )
end

function _timed_nauty(expression)
    result = @timed nauty_quotient_loop_momenta(expression)
    return result.time * 1e9, result.bytes
end

function _timed_gc(expression, workspace)
    result = @timed matrixfree_gc_quotient(expression, workspace)
    return result.time * 1e9, result.bytes
end

function _adaptive_measurement(expression, workspace)
    nauty_quotient_loop_momenta(expression)
    matrixfree_gc_quotient(expression, workspace)

    nauty_ns, nauty_memory = _timed_nauty(expression)
    gc_ns, gc_memory = _timed_gc(expression, workspace)
    ratio = gc_ns / nauty_ns

    if 0.5 <= ratio <= 2.0 || gc_memory > nauty_memory
        nauty_trial = @benchmark nauty_quotient_loop_momenta($expression) samples = 7 evals = 1
        gc_trial = @benchmark matrixfree_gc_quotient($expression, $workspace) samples = 7 evals = 1
        nauty = median(nauty_trial)
        gc = median(gc_trial)
        return nauty.time, nauty.memory, gc.time, gc.memory
    end
    return nauty_ns, nauty_memory, gc_ns, gc_memory
end

cases = corpus_cases()
@assert length(cases) == 160
@assert all(
    KC._loop_gc_capacities(expression)[2] == ACCEPTANCE_LOOP_COUNT for
    (_, expression) in cases
)
workspaces = [_workspace(expression) for (_, expression) in cases]
println("incremental-child three-loop KC acceptance corpus: ", length(cases), " cases")

semantic_failures = String[]
time_regressions = Tuple{String,Float64}[]
memory_regressions = Tuple{String,Int,Int}[]
worst_ratio = ("", 0.0)
best_ratio = ("", Inf)

@testset "incremental-child three-loop KC semantics" begin
    @test length(cases) == 160
    for (index, (label, expression)) in enumerate(cases)
        workspace = workspaces[index]
        nauty = nauty_quotient_loop_momenta(expression)
        gc = matrixfree_gc_quotient(expression, workspace)
        gc_terms = KC.loop_quotient_terms(gc)
        exact =
            gc_requotient_terms(gc) == gc_terms && gc_requotient_terms(nauty) == gc_terms
        exact || push!(semantic_failures, label)
        @test exact
    end
end

@testset "incremental-child three-loop KC performance" begin
    global worst_ratio, best_ratio
    for (index, (label, expression)) in enumerate(cases)
        workspace = workspaces[index]
        nauty_ns, nauty_memory, gc_ns, gc_memory = _adaptive_measurement(
            expression, workspace
        )
        ratio = gc_ns / nauty_ns
        ratio > worst_ratio[2] && (worst_ratio = (label, ratio))
        ratio < best_ratio[2] && (best_ratio = (label, ratio))
        gc_ns > nauty_ns && push!(time_regressions, (label, ratio))
        gc_memory > nauty_memory &&
            push!(memory_regressions, (label, gc_memory, nauty_memory))

        matrixfree_gc_quotient(expression, workspace)
        stats = _incremental_stats(workspace)
        println(
            "KC-3LOOP-ACCEPT|",
            label,
            "|gc_ns=",
            round(gc_ns; digits=1),
            "|nauty_ns=",
            round(nauty_ns; digits=1),
            "|ratio=",
            round(ratio; digits=3),
            "|gc_mem=",
            gc_memory,
            "|nauty_mem=",
            nauty_memory,
            "|vertices=",
            stats.vertices,
            "|levels=",
            stats.levels,
            "|generated=",
            stats.generated,
            "|retained=",
            stats.retained,
            "|trace_discards=",
            stats.discarded_by_trace,
            "|max_frontier=",
            stats.max_frontier,
            "|paths=",
            stats.experimental_paths,
            "|image_matches=",
            stats.image_matches,
            "|quotient_discards=",
            stats.quotient_discards,
            "|generators=",
            stats.generators,
            "|orbit_skips=",
            stats.orbit_skips,
        )
    end

    println("semantic failures: ", length(semantic_failures))
    println("confirmed time regressions: ", length(time_regressions))
    println("memory regressions: ", length(memory_regressions))
    println("best time ratio: ", round(best_ratio[2]; digits=3), "x at ", best_ratio[1])
    println("worst time ratio: ", round(worst_ratio[2]; digits=3), "x at ", worst_ratio[1])

    if !isempty(time_regressions)
        println("confirmed incremental-child/Nauty time regressions:")
        for (label, ratio) in sort(time_regressions; by=last, rev=true)
            println("  ", label, ": ", round(ratio; digits=3), "x")
        end
    end
    if !isempty(memory_regressions)
        println("incremental-child/Nauty memory regressions:")
        for (label, gc_memory, nauty_memory) in memory_regressions
            println("  ", label, ": ", gc_memory, " B vs ", nauty_memory, " B")
        end
    end

    @test isempty(semantic_failures)
    @test isempty(memory_regressions)
end
