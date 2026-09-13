"""Structural witness for one exactly reducible finite-Trotter loop frequency.

`frequency_index` is the coordinate in the canonical causal loop-frequency vector.
`loop_basis_index` is the corresponding coordinate in the source `MomentumBasis`. The latter is
retained after equal-time integration because only the loop frequency, not the later spatial
loop momentum, has been integrated.
"""
struct TrotterFrequencyWitness
    line_index::Int
    frequency_index::Int
    loop_basis_index::Int
    loop_coefficient::EnergyCoefficient
end

TrotterFrequencyWitness() = TrotterFrequencyWitness(0, 0, 0, zero(EnergyCoefficient))

"""Return whether a witness certifies an isolated equal-time loop frequency."""
function has_isolated_trotter_frequency(witness::TrotterFrequencyWitness)
    return !iszero(witness.line_index) &&
           !iszero(witness.frequency_index) &&
           !iszero(witness.loop_basis_index) &&
           !iszero(witness.loop_coefficient)
end

"""
Partial frequency state of one shifted canonical collision contribution.

The source `SpectralDispersiveTerm` and its full `MomentumBasis` are immutable provenance.
`active_frequency_indices` refers to coordinates of the canonical causal-frequency vector,
whereas `active_line_indices` records propagator factors that still participate in frequency
algebra. Equal-time integration changes neither spatial loop-momentum provenance nor the
statistical monomial.
"""
struct TrotterFrequencyState{C<:Number,S<:Statistics,E1,E2}
    source::SpectralDispersiveTerm{S,E1,E2}
    coefficient::C
    statistical::StatisticalMonomial{S}
    parameter::ParameterMonomial
    active_line_indices::Vector{Int}
    active_frequency_indices::Vector{Int}
end

function TrotterFrequencyState(
    term::ShiftedFrequencyCollisionTerm{C,S,E1,E2}
) where {C<:Number,S<:Statistics,E1,E2}
    source = frequency_source(term)
    nfrequencies = length(_loop_frequency_basis_indices(source))
    return TrotterFrequencyState{C,S,E1,E2}(
        source,
        source_coefficient(term),
        statistical_monomial(term),
        parameters(term),
        collect(eachindex(kinetic_lines(source.carrier))),
        collect(1:nfrequencies),
    )
end

frequency_source(state::TrotterFrequencyState) = state.source
source_coefficient(state::TrotterFrequencyState) = state.coefficient
statistical_monomial(state::TrotterFrequencyState) = state.statistical
parameters(state::TrotterFrequencyState) = state.parameter
momentum_basis(state::TrotterFrequencyState) = momentum_basis(state.source)
kinematic_factor(state::TrotterFrequencyState) = kinematic_factor(state.source)
active_trotter_lines(state::TrotterFrequencyState) = state.active_line_indices
active_trotter_frequencies(state::TrotterFrequencyState) = state.active_frequency_indices

function _trotter_loop_basis_index(state::TrotterFrequencyState, frequency_index::Int)
    loop_basis_indices = _loop_frequency_basis_indices(state.source)
    1 <= frequency_index <= length(loop_basis_indices) ||
        throw(BoundsError(loop_basis_indices, frequency_index))
    return loop_basis_indices[frequency_index]
end

function _trotter_frequency_incidence_count(
    state::TrotterFrequencyState, frequency_index::Int
)
    basis_index = _trotter_loop_basis_index(state, frequency_index)
    lines = kinetic_lines(state.source.carrier)
    count = 0
    for line_index in state.active_line_indices
        iszero(momentum(lines[line_index])[basis_index]) || (count += 1)
    end
    return count
end

"""
    trotter_frequency_witness(state)

Find the first shifted active line with exactly one structurally isolated active loop frequency.
A frequency is isolated when it occurs in that line and in no other active propagator factor.
If a shifted line has more than one such frequency, collapsing one propagator would leave an
unconstrained frequency volume, so no local equal-time rule is assigned.

The criterion depends only on exact affine routing and relative contour-time shift. It does not
inspect topology, perturbative order, or line count.
"""
function trotter_frequency_witness(state::TrotterFrequencyState)
    lines = kinetic_lines(state.source.carrier)

    for line_index in state.active_line_indices
        line = lines[line_index]
        iszero(regularisation_shift(line)) && continue
        routed = momentum(line)
        isolated_frequency = 0
        isolated_basis = 0
        isolated_count = 0

        for frequency_index in state.active_frequency_indices
            basis_index = _trotter_loop_basis_index(state, frequency_index)
            coefficient = routed[basis_index]
            iszero(coefficient) && continue
            _trotter_frequency_incidence_count(state, frequency_index) == 1 || continue
            isolated_count += 1
            isolated_frequency = frequency_index
            isolated_basis = basis_index
            isolated_count > 1 && break
        end

        isolated_count == 1 || continue
        coefficient = convert(EnergyCoefficient, routed[isolated_basis])
        return TrotterFrequencyWitness(
            line_index, isolated_frequency, isolated_basis, coefficient
        )
    end

    return TrotterFrequencyWitness()
end

"""Return the shifted line selected by a certified Trotter witness."""
function trotter_frequency_line(
    state::TrotterFrequencyState, witness::TrotterFrequencyWitness
)
    has_isolated_trotter_frequency(witness) ||
        throw(ArgumentError("Trotter witness does not identify an isolated frequency"))
    return kinetic_lines(state.source.carrier)[witness.line_index]
end

"""
    trotter_equal_time_factor(state, witness)

Return the exact local factor obtained by integrating the selected one-sided equal-time line.
In package conventions,

```text
A(0ˢ) = 1,
D(0ˢ) = -(i/2) sign(s),
```

multiplied by the exact routing Jacobian `1/abs(c)` for
`ω_line = c*ω_loop + ...`. A shifted statistically weighted line is deliberately not assigned a
local constant because the integral of `F(ω)A(ω)` is not fixed by the equal-time canonical
(anti)commutator alone.
"""
function trotter_equal_time_factor(
    state::TrotterFrequencyState, witness::TrotterFrequencyWitness
)::ComplexRationals
    line = trotter_frequency_line(state, witness)
    statistical_weight(line) === NoStatisticalWeight || throw(
        ArgumentError(
            "statistically weighted shifted line has no universal local equal-time factor",
        ),
    )

    shift = regularisation_shift(line)
    iszero(shift) && throw(ArgumentError("Trotter witness selected an unshifted line"))
    jacobian = inv(abs(witness.loop_coefficient))
    kind = spectral_dispersive_kinds(state.source)[witness.line_index]

    local_factor = if kind === CollisionSpectral
        one(ComplexRationals)
    elseif kind === CollisionDispersive
        convert(ComplexRationals, -(sign(Int(shift)) // 2) * im)
    else
        error("unsupported spectral/dispersive kind in Trotter reduction")
    end

    return convert(ComplexRationals, jacobian) * local_factor
end

"""
Integrate one certified isolated equal-time frequency while retaining its spatial momentum.
Only the selected propagator factor and loop-frequency coordinate are removed from active
frequency algebra; the source term and full `MomentumBasis` remain unchanged.
"""
function eliminate_trotter_frequency(
    state::TrotterFrequencyState{C,S,E1,E2}, witness::TrotterFrequencyWitness
) where {C<:Number,S<:Statistics,E1,E2}
    has_isolated_trotter_frequency(witness) ||
        throw(ArgumentError("cannot eliminate an unresolved Trotter frequency"))
    witness.line_index in state.active_line_indices ||
        throw(ArgumentError("Trotter witness line is no longer active"))
    witness.frequency_index in state.active_frequency_indices ||
        throw(ArgumentError("Trotter witness frequency is no longer active"))

    expected = trotter_frequency_witness(state)
    expected == witness || throw(
        ArgumentError(
            "Trotter witness is not the canonical isolated frequency of this state"
        ),
    )

    active_lines = Int[
        index for index in state.active_line_indices if index != witness.line_index
    ]
    active_frequencies = Int[
        index for
        index in state.active_frequency_indices if index != witness.frequency_index
    ]
    D = promote_type(C, ComplexRationals)
    coefficient =
        convert(D, state.coefficient) *
        convert(D, trotter_equal_time_factor(state, witness))
    return TrotterFrequencyState{D,S,E1,E2}(
        state.source,
        coefficient,
        state.statistical,
        state.parameter,
        active_lines,
        active_frequencies,
    )
end

"""Eliminate all consecutively exposed isolated equal-time frequencies."""
function eliminate_isolated_trotter_frequencies(state::TrotterFrequencyState)
    current = state
    while true
        witness = trotter_frequency_witness(current)
        has_isolated_trotter_frequency(witness) || return current
        current = eliminate_trotter_frequency(current, witness)
    end
end

"""Return whether an active finite-Trotter factor still requires nonlocal treatment."""
function has_unresolved_trotter_frequency(state::TrotterFrequencyState)
    lines = kinetic_lines(state.source.carrier)
    return any(
        line_index -> !iszero(regularisation_shift(lines[line_index])),
        state.active_line_indices,
    )
end

"""Return whether every loop frequency of this contribution has been integrated."""
function trotter_frequency_complete(state::TrotterFrequencyState)
    return isempty(state.active_frequency_indices)
end

"""
Build the canonical causal expression left after certified local equal-time integrations.

The denominator vectors retain the original full loop-frequency coordinate system. Coordinates
already integrated by the Trotter rule are identically absent from all surviving factors by the
isolation criterion; downstream #299/#301 reduction therefore acts only on
`active_trotter_frequencies(state)` and retains a common canonical coordinate system.
"""
function trotter_residual_causal_expression(
    state::TrotterFrequencyState{C,S}, target::FieldFamily{S}
) where {C<:Number,S<:Statistics}
    has_unresolved_trotter_frequency(state) && throw(
        ArgumentError(
            "nonisolated finite-Trotter factor remains; no ordinary causal expression exists yet",
        ),
    )

    D = promote_type(C, ComplexRationals)
    partials = CausalFrequencyTerm{D,S}[CausalFrequencyTerm(
        convert(D, state.coefficient), CausalFrequencyDenominator{S}[]
    )]
    lines = kinetic_lines(state.source.carrier)
    kinds = spectral_dispersive_kinds(state.source)

    for line_index in state.active_line_indices
        next = CausalFrequencyTerm{D,S}[]
        sizehint!(next, 2 * length(partials))
        for partial in partials
            for (factor, infinitesimal) in
                _causal_frequency_components(kinds[line_index], D)
                denominators = copy(causal_frequency_denominators(partial))
                push!(
                    denominators,
                    _causal_frequency_denominator(
                        state.source, lines[line_index], target, infinitesimal
                    ),
                )
                push!(
                    next,
                    CausalFrequencyTerm(
                        causal_frequency_coefficient(partial) * factor, denominators
                    ),
                )
            end
        end
        partials = next
    end
    return CausalFrequencyExpression(partials)
end
