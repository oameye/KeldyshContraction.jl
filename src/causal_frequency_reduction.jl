"""
Exact affine causal denominator used internally by general quasiparticle frequency reduction.

It represents `c⋅ω + E + iη0`, where `loop_coefficients` are exact rational coefficients of
the remaining loop frequencies, `energy` is the frequency-independent symbolic energy form,
and `infinitesimal` carries the exact causal prescription. The magnitude of a nonzero
infinitesimal does not affect the final Sokhotski--Plemelj limit, but retaining it exactly is
necessary because affine pole substitution can make prescriptions cancel.
"""
struct CausalFrequencyDenominator{S<:Statistics}
    loop_coefficients::Vector{EnergyCoefficient}
    energy::EnergyForm{S}
    infinitesimal::EnergyCoefficient
end

function CausalFrequencyDenominator(
    loop_coefficients::AbstractVector{<:Rational},
    energy::EnergyForm{S},
    infinitesimal::Rational,
) where {S<:Statistics}
    return CausalFrequencyDenominator{S}(
        EnergyCoefficient[convert(EnergyCoefficient, c) for c in loop_coefficients],
        energy,
        convert(EnergyCoefficient, infinitesimal),
    )
end

function Base.isequal(
    a::CausalFrequencyDenominator{S}, b::CausalFrequencyDenominator{S}
) where {S<:Statistics}
    return isequal(a.loop_coefficients, b.loop_coefficients) &&
           isequal(a.energy, b.energy) &&
           a.infinitesimal == b.infinitesimal
end
Base.:(==)(a::CausalFrequencyDenominator, b::CausalFrequencyDenominator) = isequal(a, b)
function Base.hash(denominator::CausalFrequencyDenominator, h::UInt)
    return hash(
        denominator.infinitesimal,
        hash(
            denominator.energy,
            hash(denominator.loop_coefficients, hash(CausalFrequencyDenominator, h)),
        ),
    )
end

function Base.isless(
    a::CausalFrequencyDenominator{S}, b::CausalFrequencyDenominator{S}
) where {S<:Statistics}
    length(a.loop_coefficients) == length(b.loop_coefficients) ||
        return length(a.loop_coefficients) < length(b.loop_coefficients)
    @inbounds for i in eachindex(a.loop_coefficients, b.loop_coefficients)
        a.loop_coefficients[i] == b.loop_coefficients[i] ||
            return a.loop_coefficients[i] < b.loop_coefficients[i]
    end
    isequal(a.energy, b.energy) || return isless(a.energy, b.energy)
    return a.infinitesimal < b.infinitesimal
end

"""One exact product of affine causal frequency denominators."""
struct CausalFrequencyTerm{C<:Number,S<:Statistics}
    coefficient::C
    denominators::Vector{CausalFrequencyDenominator{S}}
end

function CausalFrequencyTerm(
    coefficient::C, denominators::AbstractVector{CausalFrequencyDenominator{S}}
) where {C<:Number,S<:Statistics}
    canonical = collect(denominators)
    if !isempty(canonical)
        nfrequencies = length(first(canonical).loop_coefficients)
        basis_size = energy_basis_size(first(canonical).energy)
        for denominator in canonical
            length(denominator.loop_coefficients) == nfrequencies || throw(
                DimensionMismatch("causal denominators use different loop-frequency bases"),
            )
            energy_basis_size(denominator.energy) == basis_size ||
                throw(DimensionMismatch("causal denominators use different energy bases"))
        end
    end
    sort!(canonical)
    return CausalFrequencyTerm{C,S}(coefficient, canonical)
end

function Base.isequal(a::CausalFrequencyTerm{C,S}, b::CausalFrequencyTerm{C,S}) where {C,S}
    return a.coefficient == b.coefficient && isequal(a.denominators, b.denominators)
end
Base.:(==)(a::CausalFrequencyTerm, b::CausalFrequencyTerm) = isequal(a, b)
function Base.hash(term::CausalFrequencyTerm, h::UInt)
    return hash(term.denominators, hash(term.coefficient, hash(CausalFrequencyTerm, h)))
end

function _substitute_causal_frequency_pole(
    denominator::CausalFrequencyDenominator{S},
    pivot::CausalFrequencyDenominator{S},
    frequency_index::Int,
) where {S<:Statistics}
    pivot_coefficient = pivot.loop_coefficients[frequency_index]
    iszero(pivot_coefficient) &&
        throw(ArgumentError("causal residue pivot does not depend on selected frequency"))
    coefficient = denominator.loop_coefficients[frequency_index]
    iszero(coefficient) && return denominator

    ratio = coefficient / pivot_coefficient
    transformed = copy(denominator.loop_coefficients)
    @inbounds for i in eachindex(transformed, pivot.loop_coefficients)
        transformed[i] -= ratio * pivot.loop_coefficients[i]
    end
    transformed[frequency_index] = zero(EnergyCoefficient)
    return CausalFrequencyDenominator{S}(
        transformed,
        denominator.energy - ratio * pivot.energy,
        denominator.infinitesimal - ratio * pivot.infinitesimal,
    )
end

@inline function _causal_pole_is_upper(
    denominator::CausalFrequencyDenominator, frequency_index::Int
)
    coefficient = denominator.loop_coefficients[frequency_index]
    iszero(coefficient) &&
        throw(ArgumentError("causal pole query requires a frequency-dependent denominator"))
    iszero(denominator.infinitesimal) && throw(
        ArgumentError(
            "zero infinitesimal requires dispersive re-expansion before contour integration",
        ),
    )
    return denominator.infinitesimal / coefficient < 0
end

"""
Integrate one loop frequency of an exact product of simple causal denominators by residues.

The real-axis measure is `dω/(2π)`, and the contour is closed deterministically in the upper
half-plane. The integrand must contain at least two factors depending on the selected frequency
so the large semicircle vanishes. Zero prescriptions and coincident higher-order poles are left
for dedicated handling rather than assigned an arbitrary value.
"""
function integrate_causal_frequency(
    term::CausalFrequencyTerm{C,S}, frequency_index::Int
) where {C<:Number,S<:Statistics}
    isempty(term.denominators) &&
        throw(ArgumentError("cannot integrate a frequency-independent causal term"))
    nfrequencies = length(first(term.denominators).loop_coefficients)
    checkbounds(1:nfrequencies, frequency_index)

    dependent = Int[
        i for (i, denominator) in enumerate(term.denominators) if
        !iszero(denominator.loop_coefficients[frequency_index])
    ]
    length(dependent) >= 2 || throw(
        ArgumentError(
            "causal contour integration requires at least two frequency-dependent factors",
        ),
    )

    D = promote_type(C, ComplexRationals)
    out = CausalFrequencyTerm{D,S}[]
    for pivot_index in dependent
        pivot = term.denominators[pivot_index]
        _causal_pole_is_upper(pivot, frequency_index) || continue

        pivot_coefficient = pivot.loop_coefficients[frequency_index]
        denominators = CausalFrequencyDenominator{S}[]
        sizehint!(denominators, length(term.denominators) - 1)
        for i in eachindex(term.denominators)
            i == pivot_index && continue
            transformed = _substitute_causal_frequency_pole(
                term.denominators[i], pivot, frequency_index
            )
            if all(iszero, transformed.loop_coefficients) &&
                iszero(transformed.energy) &&
                iszero(transformed.infinitesimal)
                throw(
                    ArgumentError(
                        "coincident causal poles require a dedicated higher-order residue rule",
                    ),
                )
            end
            push!(denominators, transformed)
        end

        residue_factor = convert(D, im) / pivot_coefficient
        push!(
            out,
            CausalFrequencyTerm(
                convert(D, term.coefficient) * residue_factor, denominators
            ),
        )
    end
    return out
end

function _push_frequency_support_coefficient!(
    out::Dict{FrequencySupport{S},C}, support::FrequencySupport{S}, coefficient::C
) where {C<:Number,S<:Statistics}
    iszero(coefficient) && return out
    combined = get(out, support, zero(C)) + coefficient
    if iszero(combined)
        delete!(out, support)
    else
        out[support] = combined
    end
    return out
end

function _append_pv_branch!(
    next_states::Vector{Tuple{Vector{EnergyShell{S}},Vector{PrincipalValueSupport{S}},C}},
    shells::Vector{EnergyShell{S}},
    principal_values::Vector{PrincipalValueSupport{S}},
    coefficient::C,
    denominator::CausalFrequencyDenominator{S},
) where {C<:Number,S<:Statistics}
    pv, factor = principal_value_support(denominator.energy)
    values = copy(principal_values)
    push!(values, pv)
    push!(next_states, (copy(shells), values, coefficient * convert(C, factor)))
    return next_states
end

function _append_shell_branch!(
    next_states::Vector{Tuple{Vector{EnergyShell{S}},Vector{PrincipalValueSupport{S}},C}},
    shells::Vector{EnergyShell{S}},
    principal_values::Vector{PrincipalValueSupport{S}},
    coefficient::C,
    denominator::CausalFrequencyDenominator{S},
) where {C<:Number,S<:Statistics}
    shell, factor = energy_shell(denominator.energy)
    values = copy(shells)
    push!(values, shell)
    prescription_sign = denominator.infinitesimal > 0 ? one(C) : -one(C)
    shell_coefficient = -convert(C, im) * prescription_sign / 2
    push!(
        next_states,
        (
            values,
            copy(principal_values),
            coefficient * shell_coefficient * convert(C, factor),
        ),
    )
    return next_states
end

"""
Lower a frequency-independent causal term to canonical principal-value and shell support.

For a nonzero prescription, `1/(E+iη0) = PV(1/E) - i sign(η) A(E)/2` in the package
convention `A=i(Gᴿ-Gᴬ)`. A zero prescription is a pure principal value. Products are expanded
exactly and duplicate supports are combined.
"""
function lower_causal_frequency_term(
    term::CausalFrequencyTerm{C,S}
) where {C<:Number,S<:Statistics}
    D = promote_type(C, ComplexRationals)
    State = Tuple{Vector{EnergyShell{S}},Vector{PrincipalValueSupport{S}},D}
    states = State[(
        EnergyShell{S}[], PrincipalValueSupport{S}[], convert(D, term.coefficient)
    )]

    for denominator in term.denominators
        all(iszero, denominator.loop_coefficients) ||
            throw(ArgumentError("causal term still depends on loop frequencies"))
        iszero(denominator.energy) && throw(
            ArgumentError("zero causal energy denominator requires exceptional handling"),
        )

        next_states = State[]
        branches = iszero(denominator.infinitesimal) ? 1 : 2
        sizehint!(next_states, branches * length(states))
        for (shells, principal_values, coefficient) in states
            _append_pv_branch!(
                next_states, shells, principal_values, coefficient, denominator
            )
            if !iszero(denominator.infinitesimal)
                _append_shell_branch!(
                    next_states, shells, principal_values, coefficient, denominator
                )
            end
        end
        states = next_states
    end

    out = Dict{FrequencySupport{S},D}()
    for (shells, principal_values, coefficient) in states
        support = FrequencySupport(shells, principal_values)
        _push_frequency_support_coefficient!(out, support, coefficient)
    end
    return out
end
