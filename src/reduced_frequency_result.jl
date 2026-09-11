"""
One collision term whose strict quasiparticle shell system is linearly dependent.

`residual_support` contains the finite shell/PV/Jacobian factors obtained after exact loop-
frequency elimination while the singular dependency geometry remains in `analysis`. For a term
that also requires equal-time/Trotter regularisation this dictionary is intentionally empty and
must not be interpreted as a vanishing contribution.
"""
struct DependentFrequencyTerm{S<:Statistics,E1,E2}
    term::SpectralDispersiveTerm{S,E1,E2}
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
function requires_trotter_regularisation(term::DependentFrequencyTerm)
    return any(
        line -> !iszero(regularisation_shift(line)), kinetic_lines(term.term.carrier)
    )
end

"""
Concrete result of frequency reduction for one exact spectral/dispersive collision term.

Exactly regular quasiparticle support is stored in `regular`. Singular linearly dependent shell
systems are retained in `dependent`. Shifted equal-time structures whose evaluation requires
the Trotter prescription are retained in `trotter`. The outer result type is fixed by the input
term type and never chosen from a runtime union.

Other unsupported causal structures are deliberately *not* caught by this container: the
underlying reducer must continue to fail loudly until a mathematically derived rule exists.
"""
struct FrequencyReductionResult{S<:Statistics,E1,E2}
    regular::Dict{FrequencySupport{S},ComplexRationals}
    dependent::Vector{DependentFrequencyTerm{S,E1,E2}}
    trotter::Vector{ExceptionalFrequencyClassification{S,E1,E2}}
end

function FrequencyReductionResult(term::SpectralDispersiveTerm{S,E1,E2}) where {S,E1,E2}
    return FrequencyReductionResult{S,E1,E2}(
        Dict{FrequencySupport{S},ComplexRationals}(),
        DependentFrequencyTerm{S,E1,E2}[],
        ExceptionalFrequencyClassification{S,E1,E2}[],
    )
end

function Base.isempty(result::FrequencyReductionResult)
    return isempty(result.regular) && isempty(result.dependent) && isempty(result.trotter)
end

function _dependent_frequency_term(
    term::SpectralDispersiveTerm{S,E1,E2},
    target::FieldFamily{S},
    analysis::SpectralDependencyAnalysis{S},
) where {S<:Statistics,E1,E2}
    if any(line -> !iszero(regularisation_shift(line)), kinetic_lines(term.carrier))
        return DependentFrequencyTerm{S,E1,E2}(
            term, analysis, Dict{FrequencySupport{S},ComplexRationals}()
        )
    end
    partial = dependent_spectral_frequency_reduction(term, target)
    residual = reduce_partial_spectral_frequency(partial)
    return DependentFrequencyTerm{S,E1,E2}(term, analysis, residual)
end

"""
    reduce_frequency_term(term, target)

Reduce one collision term through every mathematically certified frequency rule while
preserving dependent-shell and equal-time/Trotter structures explicitly.

Dependent shell geometry is detected before Trotter classification because mass-shell
dependency and equal-time contour shifts are independent pieces of physics. An unshifted
dependent term is reduced as far as mathematically possible and retains the finite residual
shell/PV support multiplying its singular geometry. A term carrying both dependency and a
Trotter shift remains in `dependent` but is not assigned a residual value. Isolated
Kramers--Kronig terms reduce to the zero result. Ordinary shifted equal-time terms are retained
in `trotter`. All remaining terms are passed to the regular full-rank/causal reducer; any
unsupported causal structure continues to throw rather than being silently reclassified.
"""
function reduce_frequency_term(
    term::SpectralDispersiveTerm{S,E1,E2}, target::FieldFamily{S}
) where {S<:Statistics,E1,E2}
    result = FrequencyReductionResult(term)
    dependency = analyze_spectral_dependencies(term, target)
    if has_dependent_shell_support(dependency)
        push!(result.dependent, _dependent_frequency_term(term, target, dependency))
        return result
    end

    classification = classify_exceptional_frequency(term)
    kind = exceptional_frequency_kind(classification)
    if kind === FrequencyKramersKronigZero
        return result
    elseif kind === FrequencyTrotterRequired
        push!(result.trotter, classification)
        return result
    end

    merge!(result.regular, general_frequency_reduction(term, target))
    return result
end
