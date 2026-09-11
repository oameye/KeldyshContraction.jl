"""Kind of unresolved causal structure preserved by exact frequency reduction."""
@enum CausalExceptionalKind::UInt8 begin
    CausalCoincidentPole = 0
    CausalZeroEnergyDenominator = 1
end

"""Deterministic witness that two causal factors form one same-prescription higher-order pole."""
struct CoincidentCausalPoleWitness
    frequency_index::Int
    pivot_denominator::Int
    coincident_denominator::Int
end

CoincidentCausalPoleWitness() = CoincidentCausalPoleWitness(0, 0, 0)

function has_coincident_causal_pole(witness::CoincidentCausalPoleWitness)
    return !iszero(witness.frequency_index) &&
           !iszero(witness.pivot_denominator) &&
           !iszero(witness.coincident_denominator)
end

"""Witness for a frequency-independent causal factor whose exact energy vanishes."""
struct ZeroEnergyCausalDenominatorWitness
    denominator_index::Int
end

ZeroEnergyCausalDenominatorWitness() = ZeroEnergyCausalDenominatorWitness(0)

function has_zero_energy_causal_denominator(witness::ZeroEnergyCausalDenominatorWitness)
    return !iszero(witness.denominator_index)
end

@inline function _causal_denominator_is_zero(denominator::CausalFrequencyDenominator)
    return all(iszero, denominator.loop_coefficients) &&
           iszero(denominator.energy) &&
           iszero(denominator.infinitesimal)
end

@inline function _causal_denominator_has_zero_energy(
    denominator::CausalFrequencyDenominator
)
    return all(iszero, denominator.loop_coefficients) && iszero(denominator.energy)
end

"""
Find the first upper-half-plane causal pole whose multiplicity is greater than one.

Two factors are coincident only when eliminating the selected frequency with one as pivot makes
the other denominator identically zero in the full exact affine space, including its `i0`
prescription. This is therefore a same-pole higher-order residue problem, not an opposite-side
retarded/advanced pinch.
"""
function coincident_causal_pole_witness(term::CausalFrequencyTerm, frequency_index::Int)
    isempty(term.denominators) && return CoincidentCausalPoleWitness()
    nfrequencies = length(first(term.denominators).loop_coefficients)
    checkbounds(1:nfrequencies, frequency_index)

    dependent = Int[
        index for (index, denominator) in enumerate(term.denominators) if
        !iszero(denominator.loop_coefficients[frequency_index])
    ]
    for pivot_index in dependent
        pivot = term.denominators[pivot_index]
        _causal_pole_is_upper(pivot, frequency_index) || continue
        for denominator_index in dependent
            denominator_index == pivot_index && continue
            transformed = _substitute_causal_frequency_pole(
                term.denominators[denominator_index], pivot, frequency_index
            )
            _causal_denominator_is_zero(transformed) || continue
            return CoincidentCausalPoleWitness(
                frequency_index, pivot_index, denominator_index
            )
        end
    end
    return CoincidentCausalPoleWitness()
end

"""
Find a frequency-independent exact zero-energy causal denominator.

Such a factor is `1/(iη0)` for nonzero `η` and `PV(1/0)` for zero `η`. Neither object has a
finite strict-quasiparticle value. It must remain explicit rather than entering ordinary
Plemelj lowering. This classification is intentionally independent of whether the underlying
singularity ultimately requires a pinch/resummation interpretation.
"""
function zero_energy_causal_denominator_witness(term::CausalFrequencyTerm)
    for (index, denominator) in enumerate(term.denominators)
        _causal_denominator_has_zero_energy(denominator) || continue
        return ZeroEnergyCausalDenominatorWitness(index)
    end
    return ZeroEnergyCausalDenominatorWitness()
end

"""
Exact unresolved causal contribution exposed by frequency reduction.

`state` retains the original collision-term provenance and the active frequency/momentum split.
`dependency` records any independent/dependent spectral-shell geometry already established.
`support` contains finite shell/PV support accumulated before the exceptional causal factor.
`causal_term` is the exact remaining causal product, with all previously generated residue and
Trotter coefficients folded into its coefficient. Same-prescription higher-order poles are now
consumed by the exact derivative-residue layer; zero-energy causal denominators remain explicit.
`CausalCoincidentPole` is retained as a representation-compatible kind for previously exposed
states but is no longer emitted by the structured reducer.
"""
struct CausalExceptionalFrequencyTerm{S<:Statistics,E1,E2}
    state::FrequencyIntegrationState{S,E1,E2}
    dependency::SpectralDependencyAnalysis{S}
    support::FrequencySupport{S}
    causal_term::CausalFrequencyTerm{ComplexRationals,S}
    kind::CausalExceptionalKind
    coincident_witness::CoincidentCausalPoleWitness
    zero_energy_witness::ZeroEnergyCausalDenominatorWitness
end

function dependent_shell_support(term::CausalExceptionalFrequencyTerm)
    return dependent_shell_support(term.dependency)
end
function has_dependent_shell_support(term::CausalExceptionalFrequencyTerm)
    return has_dependent_shell_support(term.dependency)
end

function coincident_causal_pole_witness(term::CausalExceptionalFrequencyTerm)
    return term.coincident_witness
end
function zero_energy_causal_denominator_witness(term::CausalExceptionalFrequencyTerm)
    return term.zero_energy_witness
end

struct StructuredCausalReduction{S<:Statistics,E1,E2}
    regular::Dict{FrequencySupport{S},ComplexRationals}
    exceptional::Vector{CausalExceptionalFrequencyTerm{S,E1,E2}}
end

function _scaled_causal_term(
    term::CausalFrequencyTerm{C,S}, factor::ComplexRationals
) where {C<:Number,S<:Statistics}
    return CausalFrequencyTerm(
        factor * convert(ComplexRationals, term.coefficient), term.denominators
    )
end

function _push_coincident_causal_exception!(
    exceptional::Vector{CausalExceptionalFrequencyTerm{S,E1,E2}},
    state::FrequencyIntegrationState{S,E1,E2},
    dependency::SpectralDependencyAnalysis{S},
    support::FrequencySupport{S},
    term::CausalFrequencyTerm,
    factor::ComplexRationals,
    witness::CoincidentCausalPoleWitness,
) where {S<:Statistics,E1,E2}
    push!(
        exceptional,
        CausalExceptionalFrequencyTerm{S,E1,E2}(
            state,
            dependency,
            support,
            _scaled_causal_term(term, factor),
            CausalCoincidentPole,
            witness,
            ZeroEnergyCausalDenominatorWitness(),
        ),
    )
    return exceptional
end

function _push_zero_energy_causal_exception!(
    exceptional::Vector{CausalExceptionalFrequencyTerm{S,E1,E2}},
    state::FrequencyIntegrationState{S,E1,E2},
    dependency::SpectralDependencyAnalysis{S},
    support::FrequencySupport{S},
    term::CausalFrequencyTerm,
    factor::ComplexRationals,
    witness::ZeroEnergyCausalDenominatorWitness,
) where {S<:Statistics,E1,E2}
    push!(
        exceptional,
        CausalExceptionalFrequencyTerm{S,E1,E2}(
            state,
            dependency,
            support,
            _scaled_causal_term(term, factor),
            CausalZeroEnergyDenominator,
            CoincidentCausalPoleWitness(),
            witness,
        ),
    )
    return exceptional
end

"""
Reduce a partial causal problem while preserving only genuinely unresolved singular support.

Simple poles follow the frozen #299 residue/Plemelj machinery exactly. Same-prescription pole
multiplicities are consumed by the exact derivative-residue layer and then re-enter this same
recursive reduction. A frequency-independent zero-energy denominator remains singular strict-QP
support and is preserved rather than lowered to `PV(1/0)` or `1/(±i0)`.
"""
function structured_causal_frequency_reduction(
    reduction::PartialSpectralFrequencyReduction{S},
    state::FrequencyIntegrationState{S,E1,E2},
    dependency::SpectralDependencyAnalysis{S},
) where {S<:Statistics,E1,E2}
    terms = CausalFrequencyTerm{ComplexRationals,S}[
        _complex_causal_term(term) for term in reduction.causal_terms
    ]
    exceptional = CausalExceptionalFrequencyTerm{S,E1,E2}[]
    outer_factor = convert(ComplexRationals, reduction.factor) * state.factor
    nresidual = length(reduction.residual_loop_basis_indices)

    for frequency_index in 1:nresidual
        next = CausalFrequencyTerm{ComplexRationals,S}[]
        for term in terms
            zero_witness = zero_energy_causal_denominator_witness(term)
            if has_zero_energy_causal_denominator(zero_witness)
                _push_zero_energy_causal_exception!(
                    exceptional,
                    state,
                    dependency,
                    reduction.support,
                    term,
                    outer_factor,
                    zero_witness,
                )
                continue
            end

            append!(next, integrate_causal_frequency_exact(term, frequency_index))
        end
        terms = next
        isempty(terms) && break
    end

    regular = Dict{FrequencySupport{S},ComplexRationals}()
    for term in terms
        zero_witness = zero_energy_causal_denominator_witness(term)
        if has_zero_energy_causal_denominator(zero_witness)
            _push_zero_energy_causal_exception!(
                exceptional,
                state,
                dependency,
                reduction.support,
                term,
                outer_factor,
                zero_witness,
            )
            continue
        end

        lowered = lower_causal_frequency_term(term)
        for (support, coefficient) in lowered
            combined_support = _combine_frequency_support(reduction.support, support)
            _push_frequency_support_coefficient!(
                regular, combined_support, outer_factor * coefficient
            )
        end
    end
    return StructuredCausalReduction{S,E1,E2}(regular, exceptional)
end
