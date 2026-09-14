using BenchmarkTools
using KeldyshContraction
using KeldyshContraction: Classical, Quantum

const SUITE = BenchmarkGroup()

include("two_body_loss.jl")
include("two_body_scattering.jl")
include("canonicalize.jl")
include("irreducible.jl")
include("fermionic_lagrangian_sum.jl")
include("derivative_fields.jl")
include("momentum_routing.jl")
include("wigner.jl")

benchmark_two_body_loss!(SUITE)
benchmark_two_body_scattering!(SUITE)
benchmark_canonicalize!(SUITE)
benchmark_irreducible!(SUITE)
benchmark_fermionic_lagrangian_sum!(SUITE)
benchmark_derivative_fields!(SUITE)
benchmark_momentum_routing!(SUITE)
benchmark_wigner!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
