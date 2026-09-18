using Test

import KeldyshContraction as KC

include(ENV["KC_CORPUS_DEFS"])

cases = corpus_cases()
@assert length(cases) == 160

inside = Pair{String,Int}[]
outside = Pair{String,Int}[]
capacities = Int[]

for (label, expression) in cases
    graph_capacity, loop_capacity = KC._loop_gc_capacities(expression)
    @assert loop_capacity == 4
    push!(capacities, graph_capacity)
    target = graph_capacity <= 64 ? inside : outside
    push!(target, label => graph_capacity)
    println(
        "FOUR-LOOP-BOUNDARY|",
        label,
        "|vertices=",
        graph_capacity,
        "|oneword=",
        graph_capacity <= 64,
    )
end

sort!(inside; by=last)
sort!(outside; by=last)

println("four-loop cases: ", length(cases))
println("one-word eligible: ", length(inside))
println("over one-word boundary: ", length(outside))
println("minimum graph capacity: ", minimum(capacities))
println("maximum graph capacity: ", maximum(capacities))
if !isempty(inside)
    println("largest one-word case: ", last(inside).first, " => ", last(inside).second)
end
if !isempty(outside)
    println(
        "smallest over-boundary case: ", first(outside).first, " => ", first(outside).second
    )
    println("over-boundary cases:")
    for (label, capacity) in outside
        println("  ", label, ": ", capacity)
    end
end

@testset "four-loop one-word boundary inventory" begin
    @test length(inside) + length(outside) == 160
    @test all(last(case) <= 64 for case in inside)
    @test all(last(case) > 64 for case in outside)
end
