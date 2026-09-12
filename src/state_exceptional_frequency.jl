function _active_kramers_kronig_frequency_witness(state::FrequencyIntegrationState)
    term = state.term
    lines = kinetic_lines(term.carrier)
    kinds = spectral_dispersive_kinds(term)

    for basis_index in state.active_loop_basis_indices
        first_line = 0
        second_line = 0
        dependent_count = 0

        for line_index in state.active_line_indices
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
Return whether the active frequency problem contains the same isolated `A*D`
Kramers--Kronig zero certified for unreduced terms.

Only active propagator factors and active loop frequencies participate. This allows an exact
zero to become visible after an independent equal-time/Trotter loop has already been
integrated, without letting the collapsed factor remain in the proof.
"""
function has_active_kramers_kronig_zero(state::FrequencyIntegrationState)
    _, basis_index = _active_kramers_kronig_frequency_witness(state)
    return !iszero(basis_index)
end
