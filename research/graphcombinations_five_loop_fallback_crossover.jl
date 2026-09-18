using BenchmarkTools
using Test

import GraphCombinations as GC
import KeldyshContraction as KC

include(ENV["KC_CORPUS_DEFS"])
include(ENV["KC_MATRIXFREE_DEFS"])

cases = corpus_cases()
@assert length(cases) == 160

boundary_cases = [
    (label, expression, KC._loop_gc_capacities(expression)[1]) for
    (label, expression) in cases if KC._loop_gc_capacities(expression)[1] >= 60
]
@assert !isempty(boundary_cases)
@assert any(vertices == 64 for (_, _, vertices) in boundary_cases)
@assert any(vertices > 64 for (_, _, vertices) in boundary_cases)

println("five-loop GC fallback boundary cases: ", length(boundary_cases))

@testset "five-loop fallback boundary exactness and latency" begin
    for (label, expression, vertices) in boundary_cases
        workspace = _workspace(expression)
        result = matrixfree_gc_quotient(expression, workspace)
        terms = KC.loop_quotient_terms(result)
        @test gc_requotient_terms(result) == terms

        trial = @benchmark matrixfree_gc_quotient($expression, $workspace) samples = 9 evals = 1
        measured = median(trial)
        println(
            "KC-5LOOP-GC-FALLBACK|",
            label,
            "|vertices=",
            vertices,
            "|gc_ns=",
            measured.time,
            "|gc_mem=",
            measured.memory,
            "|allocs=",
            measured.allocs,
        )
    end
end
