using BenchmarkTools
using KeldyshContraction
using KeldyshContraction: Classical, Quantum

const SUITE = BenchmarkGroup()

include("two_body_loss.jl")
include("two_body_scattering.jl")

benchmark_two_body_loss!(SUITE)
benchmark_two_body_scattering!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
