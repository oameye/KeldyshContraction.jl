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
include("spectral_statistical.jl")
include("off_shell_collision.jl")
include("exceptional_frequency_classification.jl")
include("frequency_reduction.jl")

benchmark_two_body_loss!(SUITE)
benchmark_two_body_scattering!(SUITE)
benchmark_canonicalize!(SUITE)
benchmark_irreducible!(SUITE)
benchmark_fermionic_lagrangian_sum!(SUITE)
benchmark_derivative_fields!(SUITE)
benchmark_momentum_routing!(SUITE)
benchmark_wigner!(SUITE)
benchmark_spectral_statistical!(SUITE)
benchmark_off_shell_collision!(SUITE)
benchmark_exceptional_frequency_classification!(SUITE)
benchmark_frequency_reduction!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
