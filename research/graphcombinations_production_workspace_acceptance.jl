using BenchmarkTools
using Test

import KeldyshContraction as KC

include(ENV["KC_CORPUS_DEFS"])

function _timed_nauty(expression)
    result = @timed nauty_quotient_loop_momenta(expression)
    return result.time * 1e9, result.bytes
end

function _timed_prepared_gc(expression, workspace)
    result = @timed KC.quotient_loop_momenta(expression, workspace)
    return result.time * 1e9, result.bytes
end

function _adaptive_prepared_measurement(expression, workspace)
    nauty_quotient_loop_momenta(expression)
    KC.quotient_loop_momenta(expression, workspace)

    nauty_ns, nauty_memory = _timed_nauty(expression)
    gc_ns, gc_memory = _timed_prepared_gc(expression, workspace)
    ratio = gc_ns / nauty_ns

    if ratio >= 0.5 || gc_memory > nauty_memory
        nauty_trial = @benchmark nauty_quotient_loop_momenta($expression) samples = 7 evals = 1
        gc_trial = @benchmark KC.quotient_loop_momenta($expression, $workspace) samples = 7 evals = 1
        nauty = median(nauty_trial)
        gc = median(gc_trial)
        return nauty.time, nauty.memory, gc.time, gc.memory
    end
    return nauty_ns, nauty_memory, gc_ns, gc_memory
end

cases = corpus_cases()
workspaces = [KC.LoopMomentumQuotientWorkspace(expression) for (_, expression) in cases]
println("prepared production GC/Nauty corpus: ", length(cases), " generated cases")

semantic_failures = String[]
time_regressions = Tuple{String,Float64}[]
memory_regressions = Tuple{String,Int,Int}[]
worst_ratio = ("", 0.0)
best_ratio = ("", Inf)

@testset "prepared production GC semantics" begin
    for (index, (label, expression)) in enumerate(cases)
        workspace = workspaces[index]
        prepared = KC.quotient_loop_momenta(expression, workspace)
        one_shot = KC.quotient_loop_momenta(expression)
        prepared_terms = KC.loop_quotient_terms(prepared)
        exact =
            prepared_terms == KC.loop_quotient_terms(one_shot) &&
            gc_requotient_terms(prepared) == prepared_terms &&
            gc_requotient_terms(nauty_quotient_loop_momenta(expression)) == prepared_terms
        if !exact
            push!(semantic_failures, label)
        end
        @test exact
    end
end

@testset "prepared production GC performance" begin
    global worst_ratio, best_ratio
    for (index, (label, expression)) in enumerate(cases)
        workspace = workspaces[index]
        nauty_ns, nauty_memory, gc_ns, gc_memory = _adaptive_prepared_measurement(
            expression, workspace
        )
        ratio = gc_ns / nauty_ns
        if ratio > worst_ratio[2]
            worst_ratio = (label, ratio)
        end
        if ratio < best_ratio[2]
            best_ratio = (label, ratio)
        end
        if gc_ns > nauty_ns
            push!(time_regressions, (label, ratio))
        end
        if gc_memory > nauty_memory
            push!(memory_regressions, (label, gc_memory, nauty_memory))
        end
        println(
            lpad(index, 3),
            "/",
            length(cases),
            " ",
            label,
            ": time=",
            round(gc_ns / 1.0e3; digits=2),
            "/",
            round(nauty_ns / 1.0e3; digits=2),
            " μs prepared-GC/Nauty (",
            round(ratio; digits=3),
            "x), memory=",
            gc_memory,
            "/",
            nauty_memory,
            " B",
        )
    end

    println("semantic failures: ", length(semantic_failures))
    println("confirmed time regressions: ", length(time_regressions))
    println("memory regressions: ", length(memory_regressions))
    println("best time ratio: ", round(best_ratio[2]; digits=3), "x at ", best_ratio[1])
    println("worst time ratio: ", round(worst_ratio[2]; digits=3), "x at ", worst_ratio[1])

    if !isempty(time_regressions)
        println("confirmed prepared-GC time regressions:")
        for (label, ratio) in sort(time_regressions; by=last, rev=true)
            println("  ", label, ": ", round(ratio; digits=3), "x")
        end
    end
    if !isempty(memory_regressions)
        println("prepared-GC memory regressions:")
        for (label, gc_memory, nauty_memory) in memory_regressions
            println("  ", label, ": ", gc_memory, " B vs ", nauty_memory, " B")
        end
    end

    @test isempty(semantic_failures)
    @test isempty(time_regressions)
    @test isempty(memory_regressions)
end
