using BenchmarkTools
using Test

import KeldyshContraction as KC

include(ENV["KC_CORPUS_DEFS"])
include(ENV["KC_MATRIXFREE_DEFS"])

function _timed_call(f, expression)
    result = @timed f(expression)
    return result.value, result.time * 1e9, result.bytes
end

function _timed_gc(expression, workspace)
    result = @timed matrixfree_gc_quotient(expression, workspace)
    return result.value, result.time * 1e9, result.bytes
end

function _repeat_if_close(expression, workspace, nauty_ns, gc_ns)
    ratio = gc_ns / nauty_ns
    0.5 <= ratio <= 2.0 || return nauty_ns, gc_ns
    nauty_trial = @benchmark nauty_quotient_loop_momenta($expression) samples = 7 evals = 1
    gc_trial = @benchmark matrixfree_gc_quotient($expression, $workspace) samples = 7 evals = 1
    return median(nauty_trial).time, median(gc_trial).time
end

cases = corpus_cases()
@assert length(cases) == 160
overboundary = [
    (label, expression, KC._loop_gc_capacities(expression)[1]) for
    (label, expression) in cases if KC._loop_gc_capacities(expression)[1] > 64
]
@assert length(overboundary) == 4
@assert all(KC._loop_gc_capacities(expression)[2] == 5 for (_, expression, _) in overboundary)

# Compile the common code path on a cheap five-loop case before timing the four >64 cases.
warm_expression = first(cases)[2]
warm_workspace = _workspace(warm_expression)
nauty_quotient_loop_momenta(warm_expression)
matrixfree_gc_quotient(warm_expression, warm_workspace)

semantic_failures = String[]
time_ratios = Pair{String,Float64}[]
memory_ratios = Pair{String,Float64}[]

@testset "five-loop over-boundary exact fallback" begin
    @test length(overboundary) == 4
    for (label, expression, vertices) in overboundary
        workspace = _workspace(expression)

        nauty, nauty_ns, nauty_mem = _timed_call(nauty_quotient_loop_momenta, expression)
        gc, gc_ns, gc_mem = _timed_gc(expression, workspace)
        gc_terms = KC.loop_quotient_terms(gc)
        exact =
            gc_requotient_terms(gc) == gc_terms && gc_requotient_terms(nauty) == gc_terms
        exact || push!(semantic_failures, label)
        @test exact

        nauty_ns, gc_ns = _repeat_if_close(expression, workspace, nauty_ns, gc_ns)
        time_ratio = gc_ns / nauty_ns
        memory_ratio = gc_mem / nauty_mem
        push!(time_ratios, label => time_ratio)
        push!(memory_ratios, label => memory_ratio)

        println(
            "KC-5LOOP-FALLBACK|",
            label,
            "|vertices=",
            vertices,
            "|gc_ns=",
            round(gc_ns; digits=1),
            "|nauty_ns=",
            round(nauty_ns; digits=1),
            "|ratio=",
            round(time_ratio; digits=3),
            "|gc_mem=",
            gc_mem,
            "|nauty_mem=",
            nauty_mem,
            "|memory_ratio=",
            round(memory_ratio; digits=3),
        )
    end

    @test isempty(semantic_failures)
end

println("semantic failures: ", length(semantic_failures))
println("worst fallback/Nauty time ratio: ", round(maximum(last, time_ratios); digits=3))
println("worst fallback/Nauty memory ratio: ", round(maximum(last, memory_ratios); digits=3))
