"""Return whether two causal factors define the same normalized pole in one frequency."""
function _same_normalized_causal_pole(
    a::CausalFrequencyDenominator{S},
    b::CausalFrequencyDenominator{S},
    frequency_index::Int,
) where {S<:Statistics}
    coefficient_a = a.loop_coefficients[frequency_index]
    coefficient_b = b.loop_coefficients[frequency_index]
    iszero(coefficient_a) && return false
    iszero(coefficient_b) && return false

    @inbounds for index in eachindex(a.loop_coefficients, b.loop_coefficients)
        a.loop_coefficients[index] * coefficient_b ==
            b.loop_coefficients[index] * coefficient_a || return false
    end
    a.energy * coefficient_b == b.energy * coefficient_a || return false
    return a.infinitesimal * coefficient_b == b.infinitesimal * coefficient_a
end

"""
Partition all factors depending on `frequency_index` into exact normalized-pole classes.

Two denominators are in the same class iff division by their selected-frequency coefficient
would give the same complete affine denominator, including every residual loop coefficient,
`EnergyForm`, and `i0` prescription. Cross multiplication keeps the comparison exact and avoids
constructing a separate normalized symbolic representation.
"""
function causal_pole_classes(term::CausalFrequencyTerm, frequency_index::Int)
    isempty(term.denominators) && return Vector{Vector{Int}}()
    nfrequencies = length(first(term.denominators).loop_coefficients)
    checkbounds(1:nfrequencies, frequency_index)

    classes = Vector{Vector{Int}}()
    for (denominator_index, denominator) in enumerate(term.denominators)
        iszero(denominator.loop_coefficients[frequency_index]) && continue
        class_index = findfirst(classes) do pole_class
            representative = term.denominators[first(pole_class)]
            _same_normalized_causal_pole(
                denominator, representative, frequency_index
            )
        end
        if isnothing(class_index)
            push!(classes, Int[denominator_index])
        else
            push!(classes[class_index], denominator_index)
        end
    end
    return classes
end

function _push_combined_causal_term!(
    out::Vector{CausalFrequencyTerm{ComplexRationals,S}},
    coefficient::ComplexRationals,
    denominators::Vector{CausalFrequencyDenominator{S}},
) where {S<:Statistics}
    iszero(coefficient) && return out
    candidate = CausalFrequencyTerm(coefficient, denominators)
    for index in eachindex(out)
        isequal(out[index].denominators, candidate.denominators) || continue
        combined = out[index].coefficient + candidate.coefficient
        if iszero(combined)
            deleteat!(out, index)
        else
            out[index] = CausalFrequencyTerm(combined, candidate.denominators)
        end
        return out
    end
    push!(out, candidate)
    return out
end

function _differentiate_causal_product_once(
    terms::Vector{CausalFrequencyTerm{ComplexRationals,S}}, frequency_index::Int
) where {S<:Statistics}
    out = CausalFrequencyTerm{ComplexRationals,S}[]
    for term in terms
        for denominator in term.denominators
            coefficient = denominator.loop_coefficients[frequency_index]
            iszero(coefficient) && continue
            differentiated = copy(term.denominators)
            push!(differentiated, denominator)
            derivative_coefficient =
                -term.coefficient * convert(ComplexRationals, coefficient)
            _push_combined_causal_term!(out, derivative_coefficient, differentiated)
        end
    end
    return out
end

function _differentiate_causal_product(
    denominators::Vector{CausalFrequencyDenominator{S}},
    frequency_index::Int,
    order::Int,
) where {S<:Statistics}
    order >= 0 || throw(ArgumentError("causal derivative order must be nonnegative"))
    terms = CausalFrequencyTerm{ComplexRationals,S}[
        CausalFrequencyTerm(one(ComplexRationals), denominators)
    ]
    for _ in 1:order
        terms = _differentiate_causal_product_once(terms, frequency_index)
        isempty(terms) && break
    end
    return terms
end

function _higher_order_pole_residue(
    term::CausalFrequencyTerm{C,S}, pole_class::Vector{Int}, frequency_index::Int
) where {C<:Number,S<:Statistics}
    pivot = term.denominators[first(pole_class)]
    _causal_pole_is_upper(pivot, frequency_index) ||
        return CausalFrequencyTerm{ComplexRationals,S}[]

    class_mask = falses(length(term.denominators))
    coefficient_product = one(EnergyCoefficient)
    for denominator_index in pole_class
        class_mask[denominator_index] = true
        coefficient_product *=
            term.denominators[denominator_index].loop_coefficients[frequency_index]
    end

    complement = CausalFrequencyDenominator{S}[
        term.denominators[index] for index in eachindex(term.denominators) if
        !class_mask[index]
    ]
    multiplicity = length(pole_class)
    derivatives = _differentiate_causal_product(
        complement, frequency_index, multiplicity - 1
    )
    isempty(derivatives) && return CausalFrequencyTerm{ComplexRationals,S}[]

    factorial_factor = factorial(multiplicity - 1)
    residue_prefactor =
        convert(ComplexRationals, term.coefficient) * convert(ComplexRationals, im) /
        convert(ComplexRationals, factorial_factor) /
        convert(ComplexRationals, coefficient_product)

    out = CausalFrequencyTerm{ComplexRationals,S}[]
    for derivative in derivatives
        substituted = CausalFrequencyDenominator{S}[
            _substitute_causal_frequency_pole(denominator, pivot, frequency_index) for
            denominator in derivative.denominators
        ]
        _push_combined_causal_term!(
            out, residue_prefactor * derivative.coefficient, substituted
        )
    end
    return out
end

"""
Integrate one causal loop frequency with exact same-prescription pole multiplicities.

When every pole is simple this dispatches directly to the frozen simple-pole residue path.
Otherwise each distinct upper-half-plane normalized pole contributes exactly once. A class of
multiplicity `m` differentiates the complementary causal product `m-1` times and carries the
exact factor `i / ((m-1)! * prod(c_a))`. Opposite `i0` prescriptions belong to different pole
classes and are therefore never converted into a same-side higher-order residue.
"""
function integrate_causal_frequency_exact(
    term::CausalFrequencyTerm{C,S}, frequency_index::Int
) where {C<:Number,S<:Statistics}
    isempty(term.denominators) &&
        throw(ArgumentError("cannot integrate a frequency-independent causal term"))
    nfrequencies = length(first(term.denominators).loop_coefficients)
    checkbounds(1:nfrequencies, frequency_index)

    pole_classes = causal_pole_classes(term, frequency_index)
    dependent_count = sum(length, pole_classes)
    dependent_count >= 2 || throw(
        ArgumentError(
            "causal contour integration requires at least two frequency-dependent factors",
        ),
    )

    all(pole_class -> length(pole_class) == 1, pole_classes) &&
        return integrate_causal_frequency(term, frequency_index)

    out = CausalFrequencyTerm{ComplexRationals,S}[]
    for pole_class in pole_classes
        residues = _higher_order_pole_residue(term, pole_class, frequency_index)
        for residue in residues
            _push_combined_causal_term!(
                out, residue.coefficient, copy(residue.denominators)
            )
        end
    end
    return out
end
