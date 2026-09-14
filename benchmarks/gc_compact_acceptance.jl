const GC_COMPACT_SHA = "96151883a7353587197cd4375f118db7c0dd797d"

function replace_exactly_once(source::String, old::String, new::String)::String
    occurrences = length(findall(old, source))
    occurrences == 1 ||
        error("expected exactly one acceptance rewrite target, found $occurrences")
    return replace(source, old => new; count=1)
end

adapter_path = joinpath(@__DIR__, "gc_wick_adapter.jl")
adapter_source = read(adapter_path, String)
adapter_source = replace_exactly_once(
    adapter_source,
    "results, search = GC._weighted_port_matchings_with_stats(problem)",
    "results, search = GC._weighted_port_matchings_compact_with_stats(problem)",
)

println("# GraphCombinations compact acceptance SHA: ", GC_COMPACT_SHA)
println("# exact KC oracle acceptance using compact GC traversal")
include_string(Main, adapter_source, "gc_wick_adapter_compact.jl")
validate_interaction("boson_g4", L_g, 4, 9)

# Reuse the existing same-run old-vs-GC timing harness, but keep the already loaded compact
# adapter and route its isolated GC generation-stage measurement through the same compact engine.
timing_path = joinpath(@__DIR__, "gc_wick_timing.jl")
timing_source = read(timing_path, String)
timing_source = replace_exactly_once(
    timing_source, "include(\"gc_wick_adapter.jl\")\n\n", ""
)
timing_source = replace_exactly_once(
    timing_source,
    "GC._weighted_port_matchings_with_stats(first(input))[1]",
    "GC._weighted_port_matchings_compact_with_stats(first(input))[1]",
)

println("# same-run KC matcher and stage timings using compact GC traversal")
include_string(Main, timing_source, "gc_wick_timing_compact.jl")
