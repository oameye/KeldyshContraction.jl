include(joinpath(@__DIR__, "graphcombinations_propagator_gc_adapter.jl"))
include(joinpath(@__DIR__, "graphcombinations_propagator_native_adapter.jl"))
using KeldyshContraction: @qfields, Boson

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

push!(
    benchmarks,
    "2-node" => [
        KC.Contraction(ψc(KC.Out()), KC.bar(ψq)(KC.Bulk(1))),
        KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(2))),
        KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.In())),
    ],
)

push!(
    benchmarks,
    "3-node" => [
        KC.Contraction(ψc(KC.Out()), KC.bar(ψq)(KC.Bulk(2))),
        KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.Bulk(1))),
        KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.Bulk(3))),
        KC.Contraction(ψc(KC.Bulk(3)), KC.bar(ψq)(KC.Bulk(2))),
        KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(3))),
        KC.Contraction(ψc(KC.Bulk(3)), KC.bar(ψq)(KC.Bulk(1))),
        KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.In())),
    ],
)

push!(
    benchmarks,
    "4-node" => [
        KC.Contraction(ψc(KC.Out()), KC.bar(ψq)(KC.Bulk(3))),
        KC.Contraction(ψc(KC.Bulk(3)), KC.bar(ψq)(KC.Bulk(1))),
        KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(4))),
        KC.Contraction(ψc(KC.Bulk(4)), KC.bar(ψq)(KC.Bulk(2))),
        KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.Bulk(3))),
        KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.In())),
    ],
)

push!(
    benchmarks,
    "5-node" => [
        KC.Contraction(ψc(KC.Out()), KC.bar(ψq)(KC.Bulk(5))),
        KC.Contraction(ψc(KC.Bulk(5)), KC.bar(ψq)(KC.Bulk(1))),
        KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(2))),
        KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.Bulk(3))),
        KC.Contraction(ψc(KC.Bulk(3)), KC.bar(ψq)(KC.Bulk(4))),
        KC.Contraction(ψc(KC.Bulk(4)), KC.bar(ψq)(KC.Bulk(5))),
        KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(3))),
        KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.Bulk(4))),
        KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.In())),
    ],
)

push!(
    benchmarks,
    "colored-2-node" => [
        KC.Contraction(ψc(KC.Out()), KC.bar(ψq)(KC.Bulk(1))),
        KC.Contraction(ηc(KC.Bulk(1)), KC.bar(ηq)(KC.Bulk(2))),
        KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.In())),
    ],
)

push!(
    benchmarks,
    "parallel-2-node" => [
        KC.Contraction(ψc(KC.Out()), KC.bar(ψq)(KC.Bulk(1))),
        KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(2))),
        KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(2))),
        KC.Contraction(ψc(KC.Bulk(2)), KC.bar(ψq)(KC.In())),
    ],
)

push!(
    benchmarks,
    "forced-relation-fallback" => [
        KC.Contraction(ψc(KC.Bulk(1)), KC.bar(ψq)(KC.Bulk(2))),
        KC.Contraction(ηc(KC.Bulk(2)), KC.bar(ηq)(KC.Bulk(1))),
    ],
)

capacity = maximum(
    length(KC.canonicalization_positions(vs)) + length(vs) for (_, vs) in benchmarks
)
gadget_scratch = GCPhysicalCanonicalizationWorkspace(capacity)
native_scratch = GCNativePhysicalCanonicalizationWorkspace(capacity)

println()
println("production-shaped warmed canonicalization benchmark")
println("  reusable GC capacity: ", capacity)
println("  compares old subdivision gadget, native-only, native hybrid, and Nauty")
println("  times are best-of-7 batch means; allocations are warmed batch means")

for (label, vs) in benchmarks
    gadget_call = () -> gc_physical_canonicalize!(gadget_scratch, vs)
    native_call = () -> gc_native_physical_canonicalize!(native_scratch, vs)
    hybrid_call = () -> gc_native_hybrid_physical_canonicalize!(native_scratch, vs)
    nauty_call = () -> KC.canonicalize(vs)

    if label == "forced-relation-fallback"
        graph_positions = KC.canonicalization_positions(vs)
        direct_graph, vertex_colors, simple =
            gc_direct_graph!(native_scratch.direct, vs, graph_positions)
        GC.canonicalize_directed!(
            native_scratch.direct.result,
            native_scratch.direct.search,
            direct_graph,
            vertex_colors,
        )
        direct_order = GC.canonical_automorphism_order(native_scratch.direct.result)
        @assert direct_order > 1
        @assert simple
        @assert !KC.uniform_coloring(vs)

        nauty_output = nauty_call()
        @assert KC.canonicalize(gadget_call()) == nauty_output
        @assert KC.canonicalize(native_call()) == nauty_output
        @assert KC.canonicalize(hybrid_call()) == nauty_output
        println(
            "FALLBACK\tforced-relation-fallback\tdirect_order=",
            direct_order,
            "\tsimple=",
            simple,
            "\tuniform=false\trelation_fallback=true",
        )
    end

    for _ in 1:100
        gadget_call()
        native_call()
        hybrid_call()
        nauty_call()
    end

    repetitions = label == "5-node" ? 200 : 1000
    gadget_time = best_batch_time(gadget_call, repetitions)
    native_time = best_batch_time(native_call, repetitions)
    hybrid_time = best_batch_time(hybrid_call, repetitions)
    nauty_time = best_batch_time(nauty_call, repetitions)
    gadget_alloc = allocation_per_call(gadget_call, repetitions)
    native_alloc = allocation_per_call(native_call, repetitions)
    hybrid_alloc = allocation_per_call(hybrid_call, repetitions)
    nauty_alloc = allocation_per_call(nauty_call, repetitions)

    println(
        "MICRO\t",
        label,
        "\tgadget_over_nauty=",
        round(gadget_time / nauty_time; digits=3),
        "\tnative_over_nauty=",
        round(native_time / nauty_time; digits=3),
        "\thybrid_over_nauty=",
        round(hybrid_time / nauty_time; digits=3),
        "\tnative_over_gadget=",
        round(native_time / gadget_time; digits=3),
        "\thybrid_over_gadget=",
        round(hybrid_time / gadget_time; digits=3),
        "\tgadget_us=",
        round(gadget_time * 1e6; digits=3),
        "\tnative_us=",
        round(native_time * 1e6; digits=3),
        "\thybrid_us=",
        round(hybrid_time * 1e6; digits=3),
        "\tnauty_us=",
        round(nauty_time * 1e6; digits=3),
        "\tgadget_B=",
        round(gadget_alloc; digits=1),
        "\tnative_B=",
        round(native_alloc; digits=1),
        "\thybrid_B=",
        round(hybrid_alloc; digits=1),
        "\tnauty_B=",
        round(nauty_alloc; digits=1),
    )
end
