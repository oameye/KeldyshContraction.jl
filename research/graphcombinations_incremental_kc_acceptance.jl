using BenchmarkTools
using Test

import GraphCombinations as GC
import KeldyshContraction as KC

include(ENV["KC_CORPUS_DEFS"])
include(ENV["KC_MATRIXFREE_DEFS"])

const ACCEPTANCE_LOOP_COUNT = 2

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

function _measure_pair_15(expression, workspace)
    nauty_quotient_loop_momenta(expression)
    matrixfree_gc_quotient(expression, workspace)
    nauty_trial = @benchmark nauty_quotient_loop_momenta($expression) samples = 15 evals = 1
    gc_trial = @benchmark matrixfree_gc_quotient($expression, $workspace) samples = 15 evals =
        1
    return median(nauty_trial), median(gc_trial)
end

function _measure_pair_41(expression, workspace)
    nauty_trial = @benchmark nauty_quotient_loop_momenta($expression) samples = 41 evals = 1
    gc_trial = @benchmark matrixfree_gc_quotient($expression, $workspace) samples = 41 evals =
        1
    return median(nauty_trial), median(gc_trial)
end

function _confirmed_measurement(expression, workspace)
    nauty, gc = _measure_pair_15(expression, workspace)
    if gc.time > nauty.time || gc.memory > nauty.memory
        nauty, gc = _measure_pair_41(expression, workspace)
    end
    return nauty, gc
end

cases = corpus_cases()
workspaces = [_workspace(expression) for (_, expression) in cases]
println("incremental-child KC acceptance corpus: ", length(cases), " cases")

semantic_failures = String[]
time_regressions = Tuple{String,Float64}[]
memory_regressions = Tuple{String,Int,Int}[]
worst_ratio = ("", 0.0)
best_ratio = ("", Inf)

@testset "incremental-child KC semantics" begin
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

@testset "incremental-child KC performance" begin
    global worst_ratio, best_ratio
    for (index, (label, expression)) in enumerate(cases)
        workspace = workspaces[index]
        nauty, gc = _confirmed_measurement(expression, workspace)
        ratio = gc.time / nauty.time
        ratio > worst_ratio[2] && (worst_ratio = (label, ratio))
        ratio < best_ratio[2] && (best_ratio = (label, ratio))
        gc.time > nauty.time && push!(time_regressions, (label, ratio))
        gc.memory > nauty.memory &&
            push!(memory_regressions, (label, gc.memory, nauty.memory))

        matrixfree_gc_quotient(expression, workspace)
        stats = _incremental_stats(workspace)
        println(
            "KC-ACCEPT|",
            label,
            "|gc_ns=",
            gc.time,
            "|nauty_ns=",
            nauty.time,
            "|ratio=",
            round(ratio; digits=3),
            "|gc_mem=",
            gc.memory,
            "|nauty_mem=",
            nauty.memory,
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
