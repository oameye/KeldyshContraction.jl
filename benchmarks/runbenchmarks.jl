using BenchmarkTools
using KeldyshContraction
using KeldyshContraction: Classical, Quantum
using KeldyshContraction:
    KineticLineKind,
    KineticRetarded,
    KineticAdvanced,
    KineticSpectral,
    StatisticalWeight,
    NoStatisticalWeight,
    DistributionWeight,
    KineticLine,
    KineticMonomial,
    KineticTerm,
    KineticExpression,
    kinetic_lines,
    kinetic_line_kind,
    statistical_weight,
    regularisation_shift,
    SpectralDispersiveKind,
    CollisionDispersive,
    CollisionSpectral,
    SpectralDispersiveTerm,
    SpectralDispersiveExpression,
    spectral_dispersive_kinds,
    DispersionAtom,
    EnergyForm,
    energy_terms,
    energy_basis_size,
    energy_shell,
    principal_value_support,
    FullRankFrequencyReduction,
    loop_frequency_count,
    spectral_frequency_rank,
    full_rank_frequency_reduction,
    loop_frequency_basis_indices,
    loop_frequency_energies,
    spectral_pivot_lines,
    ExceptionalFrequencyKind,
    FrequencyUnresolved,
    FrequencyKramersKronigZero,
    FrequencyTrotterRequired,
    ExceptionalFrequencyClassification,
    classify_exceptional_frequency,
    exceptional_frequency_kind,
    exceptional_frequency_term,
    exceptional_frequency_witness_lines,
    exceptional_frequency_loop_basis_index,
    StatisticalAtom,
    StatisticalMonomial,
    StatisticalPolynomial,
    statistical_family,
    OccupationAtom,
    OccupationMonomial,
    OccupationPolynomial,
    occupation_family,
    occupation_statistics_sign,
    occupation_collision_factor,
    occupation_substitute,
    occupation_collision_polynomial,
    ReducedCollisionSector,
    ReducedDependentCollisionSector,
    ReducedCausalCollisionSector,
    ReducedTrotterCollisionTerm,
    active_loop_basis_indices,
    causal_exceptional_kind,
    causal_denominators,
    source_coefficient,
    statistical_monomial,
    frequency_branch,
    LoopMomentumTransform,
    loop_transform_matrix,
    external_momentum_index,
    loop_permutation_transform,
    transform_loop_momenta

const SUITE = BenchmarkGroup()

include("two_body_loss.jl")
include("two_body_scattering.jl")
include("canonicalize.jl")
include("irreducible.jl")
include("fermionic_lagrangian_sum.jl")
include("fermionic_pwave_loss.jl")
include("derivative_fields.jl")
include("momentum_routing.jl")
include("wigner.jl")
include("spectral_statistical.jl")
include("off_shell_collision.jl")
include("exceptional_frequency_classification.jl")
include("frequency_reduction.jl")
include("collision_reduction.jl")

benchmark_two_body_loss!(SUITE)
benchmark_two_body_scattering!(SUITE)
benchmark_canonicalize!(SUITE)
benchmark_irreducible!(SUITE)
benchmark_fermionic_lagrangian_sum!(SUITE)
benchmark_fermionic_pwave_loss!(SUITE)
benchmark_derivative_fields!(SUITE)
benchmark_momentum_routing!(SUITE)
benchmark_wigner!(SUITE)
benchmark_spectral_statistical!(SUITE)
benchmark_off_shell_collision!(SUITE)
benchmark_exceptional_frequency_classification!(SUITE)
benchmark_frequency_reduction!(SUITE)
benchmark_collision_reduction!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
