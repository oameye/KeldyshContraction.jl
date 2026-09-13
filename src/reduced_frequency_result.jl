"""
One collision term whose strict quasiparticle shell system is linearly dependent.

`state` records the exact frequency-integration state at which the dependency was established.
`residual_support` contains every mathematically finite shell/PV/Jacobian factor multiplying
the singular geometry. A remaining shifted active factor keeps this dictionary empty until its
Trotter prescription is resolved; a shift already eliminated from `state` is not treated as an
outstanding regularisation requirement.
"""
struct DependentFrequencyTerm{S<:Statistics,E1,E2}
    state::FrequencyIntegrationState{S,E1,E2}
    analysis::SpectralDependencyAnalysis{S}
    residual_support::Dict{FrequencySupport{S},ComplexRationals}
end

function dependent_shell_support(term::DependentFrequencyTerm)
    return dependent_shell_support(term.analysis)
end
function constraint_scale_factor(term::DependentFrequencyTerm)
    return constraint_scale_factor(term.analysis)
end
dependent_residual_support(term::DependentFrequencyTerm) = term.residual_support

function _has_active_trotter_shift(state::FrequencyIntegrationState)
    lines = kinetic_lines(state.term.carrier)
    return any(
        line_index -> !iszero(regularisation_shift(lines[line_index])),
        state.active_line_indices,
    )
end

function requires_trotter_regularisation(term::DependentFrequencyTerm)
    return _has_active_trotter_shift(term.state)
end

"""
Concrete result of frequency reduction for one exact spectral/dispersive collision term.

Exactly regular quasiparticle support is stored in `regular`. Singular linearly dependent shell
systems with a finite residual quotient are retained in `dependent`. Same-prescription
higher-order causal poles and frequency-independent zero-energy causal denominators are retained
in `causal`. A `trotter` entry stores the partially integrated frequency state only when shifted
active structure remains after every certified isolated equal-time loop has been eliminated.
The outer result type is fixed by the input term type and never selected through a runtime union.

The four branches are physical categories, not necessarily mutually exclusive decompositions:
a dispersive causal expansion may produce both finite regular support and an unresolved causal
branch. No unresolved causal branch is assigned a numerical value.
"""
struct FrequencyReductionResult{S<:Statistics,E1,E2}
    regular::Dict{FrequencySupport{S},ComplexRationals}
    dependent::Vector{DependentFrequencyTerm{S,E1,E2}}
    causal::Vector{CausalExceptionalFrequencyTerm{S,E1,E2}}
    trotter::Vector{FrequencyIntegrationState{S,E1,E2}}
end

function FrequencyReductionResult(term::SpectralDispersiveTerm{S,E1,E2}) where {S,E1,E2}
    return FrequencyReductionResult{S,E1,E2}(
        Dict{FrequencySupport{S},ComplexRationals}(),
        DependentFrequencyTerm{S,E1,E2}[],
        CausalExceptionalFrequencyTerm{S,E1,E2}[],
        FrequencyIntegrationState{S,E1,E2}[],
    )
end

function Base.isempty(result::FrequencyReductionResult)
    return isempty(result.regular) &&
           isempty(result.dependent) &&
           isempty(result.causal) &&
           isempty(result.trotter)
end

function _reduce_dependent_state!(
    result::FrequencyReductionResult{S,E1,E2},
    state::FrequencyIntegrationState{S,E1,E2},
    target::FieldFamily{S},
    dependency::SpectralDependencyAnalysis{S},
) where {S<:Statistics,E1,E2}
    partial = _active_dependent_spectral_frequency_reduction(state, target)
    structured = structured_causal_frequency_reduction(partial, state, dependency)
    if !isempty(structured.regular)
        push!(
            result.dependent,
            DependentFrequencyTerm{S,E1,E2}(state, dependency, structured.regular),
        )
    end
    append!(result.causal, structured.exceptional)
    return result
end

function _reduce_regular_state!(
    result::FrequencyReductionResult{S,E1,E2},
    state::FrequencyIntegrationState{S,E1,E2},
    target::FieldFamily{S},
    dependency::SpectralDependencyAnalysis{S},
) where {S<:Statistics,E1,E2}
    if spectral_frequency_rank(state) == loop_frequency_count(state)
        merge!(result.regular, general_frequency_reduction(state, target))
        return result
    end

    partial = _active_partial_spectral_frequency_reduction(state, target)
    structured = structured_causal_frequency_reduction(partial, state, dependency)
    merge!(result.regular, structured.regular)
    append!(result.causal, structured.exceptional)
    return result
end

function _reduce_unshifted_frequency_term!(
    result::FrequencyReductionResult{S,E1,E2},
    term::SpectralDispersiveTerm{S,E1,E2},
    target::FieldFamily{S},
) where {S<:Statistics,E1,E2}
    state = FrequencyIntegrationState(term)
    dependency = analyze_spectral_dependencies(state, target)
    if has_dependent_shell_support(dependency)
        return _reduce_dependent_state!(result, state, target, dependency)
    end

    classification = classify_exceptional_frequency(term)
    kind = exceptional_frequency_kind(classification)
    if kind === FrequencyKramersKronigZero
        return result
    elseif kind === FrequencyTrotterRequired
        error("unshifted frequency path received a Trotter-required term")
    end

    return _reduce_regular_state!(result, state, target, dependency)
end

function _reduce_shifted_frequency_term!(
    result::FrequencyReductionResult{S,E1,E2},
    term::SpectralDispersiveTerm{S,E1,E2},
    target::FieldFamily{S},
) where {S<:Statistics,E1,E2}
    state = eliminate_isolated_trotter_frequencies(FrequencyIntegrationState(term))
    if _has_active_trotter_shift(state)
        push!(result.trotter, state)
        return result
    end

    dependency = analyze_spectral_dependencies(state, target)
    if has_dependent_shell_support(dependency)
        return _reduce_dependent_state!(result, state, target, dependency)
    end

    has_active_kramers_kronig_zero(state) && return result
    return _reduce_regular_state!(result, state, target, dependency)
end

"""
    reduce_frequency_term(term, target)

Reduce one collision term through every mathematically certified frequency rule while keeping
ordinary shell/PV reduction, dependent-shell geometry, equal-time contour physics, and
unresolved causal singularities explicit.

A shifted term first eliminates every structurally isolated equal-time loop using the exact
one-sided Trotter rule. Those frequency integrations remove only active frequency variables and
propagator factors; their spatial loop momenta remain in the original momentum basis. Dependency
analysis is then repeated on the active shell system, so a spectral factor already collapsed at
equal time cannot spuriously participate in a pinch relation.

After Trotter elimination, shifted and originally unshifted terms share the same active-state
frequency reducer. Full-rank regular support retains the exact fast path. Every rank-deficient
causal problem goes through `structured_causal_frequency_reduction`, so same-prescription
higher-order poles and frequency-independent zero-energy denominators are preserved identically
regardless of whether a Trotter factor was present upstream. Isolated Kramers--Kronig zeros are
proved before regular causal reduction. No Trotter shift is ever used as a mass-shell width.
"""
function reduce_frequency_term(
    term::SpectralDispersiveTerm{S,E1,E2}, target::FieldFamily{S}
) where {S<:Statistics,E1,E2}
    result = FrequencyReductionResult(term)
    shifted = any(line -> !iszero(regularisation_shift(line)), kinetic_lines(term.carrier))
    if shifted
        return _reduce_shifted_frequency_term!(result, term, target)
    end
    return _reduce_unshifted_frequency_term!(result, term, target)
end
