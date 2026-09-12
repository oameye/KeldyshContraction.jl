"""Classification of frequency structures not handled by the generic full-rank reducer."""
@enum ExceptionalFrequencyKind::UInt8 begin
    FrequencyUnresolved = 0
    FrequencyKramersKronigZero = 1
    FrequencyTrotterRequired = 2
end

"""
Concrete structural classification for a frequency-reduction exception.

`term` retains the exact routed spectral/dispersive contribution. `spectral_rank` and
`loop_count` record why the generic reduction boundary was reached. For
`FrequencyKramersKronigZero`, `witness_lines == (i, j)` stores the spectral and dispersive
line indices, respectively, and `loop_basis_index` identifies the isolated loop-frequency
variable whose integration proves the cancellation. For `FrequencyTrotterRequired`, the
first witness index identifies the first shifted equal-time line. Zero witness indices mean
that no stronger structural proof is currently available.
"""
struct ExceptionalFrequencyClassification{S<:Statistics,E1,E2}
    term::SpectralDispersiveTerm{S,E1,E2}
    kind::ExceptionalFrequencyKind
    spectral_rank::Int
    loop_count::Int
    witness_lines::NTuple{2,Int}
    loop_basis_index::Int
end

"""Return the exceptional frequency classification tag."""
function exceptional_frequency_kind(classification::ExceptionalFrequencyClassification)
    return classification.kind
end

"""Return the original exact spectral/dispersive term retained by the classification."""
function exceptional_frequency_term(classification::ExceptionalFrequencyClassification)
    return classification.term
end

"""
Return the structural witness line indices.

For `FrequencyKramersKronigZero` the tuple is `(spectral_line, dispersive_line)`. For
`FrequencyTrotterRequired` the first entry identifies a shifted equal-time line and the second
entry is zero. `(0, 0)` denotes an unresolved exception.
"""
function exceptional_frequency_witness_lines(
    classification::ExceptionalFrequencyClassification
)
    return classification.witness_lines
end

"""Return the loop-basis index used by a Kramers--Kronig zero proof, or zero otherwise."""
function exceptional_frequency_loop_basis_index(
    classification::ExceptionalFrequencyClassification
)
    return classification.loop_basis_index
end

function Base.isequal(
    a::ExceptionalFrequencyClassification{S,E1,E2},
    b::ExceptionalFrequencyClassification{S,E1,E2},
) where {S<:Statistics,E1,E2}
    return isequal(a.term, b.term) &&
           a.kind === b.kind &&
           a.spectral_rank == b.spectral_rank &&
           a.loop_count == b.loop_count &&
           a.witness_lines == b.witness_lines &&
           a.loop_basis_index == b.loop_basis_index
end
function Base.:(==)(
    a::ExceptionalFrequencyClassification, b::ExceptionalFrequencyClassification
)
    return isequal(a, b)
end
function Base.hash(classification::ExceptionalFrequencyClassification, h::UInt)
    return hash(
        (
            classification.term,
            classification.kind,
            classification.spectral_rank,
            classification.loop_count,
            classification.witness_lines,
            classification.loop_basis_index,
        ),
        hash(ExceptionalFrequencyClassification, h),
    )
end

@inline function _same_causal_frequency_denominator(
    a::KineticLine{S}, b::KineticLine{S}
) where {S<:Statistics}
    # After homogeneous Fourier/Wigner lowering, the causal frequency denominator is fixed by
    # the physical species, routed momentum, and equal-time shift. Coordinate-graph endpoints
    # are provenance only: oppositely oriented edges can carry the same frequency denominator.
    return isequal(a.family, b.family) &&
           a.regularisation_shift == b.regularisation_shift &&
           isequal(momentum(a), momentum(b))
end

function _first_shifted_frequency_line(term::SpectralDispersiveTerm)
    for (line_index, line) in enumerate(kinetic_lines(term.carrier))
        iszero(regularisation_shift(line)) || return line_index
    end
    return 0
end

function _kramers_kronig_frequency_witness(term::SpectralDispersiveTerm)
    lines = kinetic_lines(term.carrier)
    kinds = spectral_dispersive_kinds(term)

    for basis_index in _loop_frequency_basis_indices(term)
        first_line = 0
        second_line = 0
        dependent_count = 0

        for line_index in eachindex(lines)
            iszero(momentum(lines[line_index])[basis_index]) && continue
            dependent_count += 1
            if dependent_count == 1
                first_line = line_index
            elseif dependent_count == 2
                second_line = line_index
            else
                break
            end
        end

        dependent_count == 2 || continue
        first_kind = kinds[first_line]
        second_kind = kinds[second_line]
        first_kind === second_kind && continue

        spectral_line = first_kind === CollisionSpectral ? first_line : second_line
        dispersive_line = first_kind === CollisionDispersive ? first_line : second_line
        statistical_weight(lines[dispersive_line]) === NoStatisticalWeight || continue
        _same_causal_frequency_denominator(lines[spectral_line], lines[dispersive_line]) ||
            continue

        return (spectral_line, dispersive_line), basis_index
    end

    return (0, 0), 0
end

"""
    classify_exceptional_frequency(term)

Classify frequency structures that require a rule beyond `full_rank_frequency_reduction`.
The result type is fixed by the input term type; the runtime outcome is stored in an
`ExceptionalFrequencyKind` tag.

The current Kramers--Kronig zero proof is deliberately conservative. It is accepted only when
one loop-frequency variable occurs in exactly two factors, those factors are one spectral and
one dispersive copy of the same routed causal denominator, and no equal-time shift is present.
The denominator identity is defined in the post-Wigner kinetic variables (field family, routed
momentum, and shift), not by the coordinate-graph endpoint labels retained as provenance. All
other factors are therefore independent of the witnessed frequency and factor out of the
frequency integral. Within the quasiparticle statistical convention, any distribution weight
on the spectral factor is likewise frequency-independent and factors out, leaving the causal
orthogonality integral `∫ A(ω) D(ω) dω = 0`.

Any nonzero equal-time regularisation shift is classified as `FrequencyTrotterRequired` before
Kramers--Kronig analysis. Rank deficiency, unmatched spectral/dispersive products, or frequency
integrals involving additional dependent factors remain `FrequencyUnresolved`; they are never
silently converted to zero or principal-value support.
"""
function classify_exceptional_frequency(
    term::SpectralDispersiveTerm{S,E1,E2}
) where {S<:Statistics,E1,E2}
    rank = spectral_frequency_rank(term)
    nloops = loop_frequency_count(term)

    shifted_line = _first_shifted_frequency_line(term)
    if !iszero(shifted_line)
        return ExceptionalFrequencyClassification{S,E1,E2}(
            term, FrequencyTrotterRequired, rank, nloops, (shifted_line, 0), 0
        )
    end

    witness_lines, loop_basis_index = _kramers_kronig_frequency_witness(term)
    if !iszero(loop_basis_index)
        return ExceptionalFrequencyClassification{S,E1,E2}(
            term, FrequencyKramersKronigZero, rank, nloops, witness_lines, loop_basis_index
        )
    end

    return ExceptionalFrequencyClassification{S,E1,E2}(
        term, FrequencyUnresolved, rank, nloops, (0, 0), 0
    )
end
