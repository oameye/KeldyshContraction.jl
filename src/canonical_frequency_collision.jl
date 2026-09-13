"""
Frequency-independent physical identity of one canonical collision-frequency expression.

Topology is deliberately absent. Terms can combine only when parameter provenance, exact momentum
routing, derivative kinematics, and the complete statistical monomial agree. This makes exact
cancellations visible before any contour decision while keeping the later spatial loop-momentum
quotient as a separate operation.
"""
struct CanonicalFrequencySector{S<:Statistics}
    parameter::ParameterMonomial
    basis::MomentumBasis
    external_momentum::MomentumVariable
    kinematic::MomentumPolynomial{ComplexRationals}
    statistical::StatisticalMonomial{S}
end

statistics(::CanonicalFrequencySector{S}) where {S<:Statistics} = S
parameters(sector::CanonicalFrequencySector) = sector.parameter
momentum_basis(sector::CanonicalFrequencySector) = sector.basis
external_wigner_momentum(sector::CanonicalFrequencySector) = sector.external_momentum
kinematic_factor(sector::CanonicalFrequencySector) = sector.kinematic
statistical_monomial(sector::CanonicalFrequencySector) = sector.statistical

function Base.isequal(
    a::CanonicalFrequencySector{S}, b::CanonicalFrequencySector{S}
) where {S}
    return isequal(a.parameter, b.parameter) &&
           isequal(a.basis, b.basis) &&
           isequal(a.external_momentum, b.external_momentum) &&
           isequal(a.kinematic, b.kinematic) &&
           isequal(a.statistical, b.statistical)
end
Base.:(==)(a::CanonicalFrequencySector, b::CanonicalFrequencySector) = isequal(a, b)
function Base.hash(sector::CanonicalFrequencySector, h::UInt)
    h = hash(CanonicalFrequencySector, h)
    h = hash(sector.parameter, h)
    h = hash(sector.basis, h)
    h = hash(sector.external_momentum, h)
    h = hash(sector.kinematic, h)
    return hash(sector.statistical, h)
end

"""
One collision contribution carrying a finite contour-time shift.

No frequency reduction is attempted here. The exact source term, numerical coefficient,
statistical monomial, and parameter provenance are retained for the later regulator layer.
"""
struct ShiftedFrequencyCollisionTerm{C<:Number,S<:Statistics,E1,E2}
    source::SpectralDispersiveTerm{S,E1,E2}
    coefficient::C
    statistical::StatisticalMonomial{S}
    parameter::ParameterMonomial
end

source_coefficient(term::ShiftedFrequencyCollisionTerm) = term.coefficient
statistical_monomial(term::ShiftedFrequencyCollisionTerm) = term.statistical
parameters(term::ShiftedFrequencyCollisionTerm) = term.parameter
frequency_source(term::ShiftedFrequencyCollisionTerm) = term.source

"""
Complete collision expression in canonical affine causal frequency form.

`expressions` contains only unshifted sectors. All exact statistical factors, including the
external `F_target(k)` multiplying the Kadanoff--Baym distribution coefficient, are assembled
before causal frequency expressions are merged. `shifted` retains finite-Trotter contributions
without assigning them a contour value.
"""
struct CanonicalFrequencyCollision{
    C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext
}
    expressions::Dict{CanonicalFrequencySector{S},CausalFrequencyExpression{C,S}}
    shifted::Vector{ShiftedFrequencyCollisionTerm{C,S,E1,E2}}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::CanonicalFrequencyCollision{C,S,O}) where {C,S,O} = O
statistics(::CanonicalFrequencyCollision{C,S}) where {C,S} = S
target_family(collision::CanonicalFrequencyCollision) = collision.target
parameters(collision::CanonicalFrequencyCollision) = collision.parameter
gradient_order(::CanonicalFrequencyCollision{C,S,O,E1,E2,G}) where {C,S,O,E1,E2,G} = Val(G)
wigner_context(collision::CanonicalFrequencyCollision) = collision.context
function canonical_frequency_expressions(collision::CanonicalFrequencyCollision)
    return collision.expressions
end
shifted_frequency_terms(collision::CanonicalFrequencyCollision) = collision.shifted
function Base.isempty(collision::CanonicalFrequencyCollision)
    return isempty(collision.expressions) && isempty(collision.shifted)
end

@inline function _canonical_frequency_coefficient_type(::Type{C}) where {C<:Number}
    return promote_type(C, ComplexRationals)
end

function _internal_statistical_monomial(
    term::SpectralDispersiveTerm{S}
) where {S<:Statistics}
    atoms = StatisticalAtom{S}[]
    for line in kinetic_lines(term.carrier)
        statistical_weight(line) === DistributionWeight || continue
        push!(atoms, StatisticalAtom{S}(line.family, momentum(line)))
    end
    return StatisticalMonomial(atoms)
end

function _external_statistical_atom(
    term::SpectralDispersiveTerm{S}, target::FieldFamily{S}
) where {S<:Statistics}
    basis = momentum_basis(term)
    external_index = _external_frequency_basis_index(term)
    return StatisticalAtom{S}(target, basis_momentum(basis, external_index))
end

function _collision_statistical_monomial(
    term::SpectralDispersiveTerm{S}, target::FieldFamily{S}, include_external::Bool
) where {S<:Statistics}
    internal = _internal_statistical_monomial(term)
    include_external || return internal
    external = StatisticalMonomial(
        StatisticalAtom{S}[_external_statistical_atom(term, target)]
    )
    return internal * external
end

function _canonical_frequency_sector(
    term::SpectralDispersiveTerm{S},
    statistical::StatisticalMonomial{S},
    parameter::ParameterMonomial,
) where {S<:Statistics}
    return CanonicalFrequencySector{S}(
        parameter,
        momentum_basis(term),
        external_wigner_momentum(term),
        kinematic_factor(term),
        statistical,
    )
end

function _causal_frequency_denominator(
    term::SpectralDispersiveTerm{S},
    line::KineticLine{S},
    target::FieldFamily{S},
    infinitesimal::EnergyCoefficient,
) where {S<:Statistics}
    basis = momentum_basis(term)
    external_index = _external_frequency_basis_index(term)
    loop_indices = _loop_frequency_basis_indices(term)
    routed = momentum(line)
    target_momentum = basis_momentum(basis, external_index)
    target_energy = EnergyForm(DispersionAtom(target, target_momentum))
    energy =
        routed[external_index] * target_energy -
        EnergyForm(DispersionAtom(line.family, routed))
    coefficients = EnergyCoefficient[routed[index] for index in loop_indices]
    return CausalFrequencyDenominator(coefficients, energy, infinitesimal)
end

@inline function _causal_frequency_components(
    kind::SpectralDispersiveKind, ::Type{D}
) where {D<:Number}
    if kind === CollisionDispersive
        half = convert(D, 1 // 2)
        return ((half, one(EnergyCoefficient)), (half, -one(EnergyCoefficient)))
    end
    kind === CollisionSpectral || error("unsupported spectral/dispersive collision kind")
    return (
        (convert(D, im), one(EnergyCoefficient)), (-convert(D, im), -one(EnergyCoefficient))
    )
end

function _causal_frequency_expression(
    term::SpectralDispersiveTerm{S,E1,E2}, source_coefficient::C, target::FieldFamily{S}
) where {C<:Number,S<:Statistics,E1,E2}
    D = _canonical_frequency_coefficient_type(C)
    partials = CausalFrequencyTerm{D,S}[CausalFrequencyTerm(
        convert(D, source_coefficient), CausalFrequencyDenominator{S}[]
    )]
    lines = kinetic_lines(term.carrier)
    kinds = spectral_dispersive_kinds(term)

    for line_index in eachindex(lines)
        next = CausalFrequencyTerm{D,S}[]
        sizehint!(next, 2 * length(partials))
        for partial in partials
            for (factor, infinitesimal) in
                _causal_frequency_components(kinds[line_index], D)
                denominators = copy(causal_frequency_denominators(partial))
                push!(
                    denominators,
                    _causal_frequency_denominator(
                        term, lines[line_index], target, infinitesimal
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

function _has_trotter_shift(term::SpectralDispersiveTerm)
    return any(line -> !iszero(regularisation_shift(line)), kinetic_lines(term.carrier))
end

function _merge_causal_frequency_expression!(
    expressions::Dict{CanonicalFrequencySector{S},CausalFrequencyExpression{D,S}},
    sector::CanonicalFrequencySector{S},
    contribution::CausalFrequencyExpression{D,S},
) where {D<:Number,S<:Statistics}
    isempty(contribution) && return expressions
    if haskey(expressions, sector)
        combined_terms = vcat(
            causal_frequency_terms(expressions[sector]),
            causal_frequency_terms(contribution),
        )
        combined = CausalFrequencyExpression(combined_terms)
        if isempty(combined)
            delete!(expressions, sector)
        else
            expressions[sector] = combined
        end
    else
        expressions[sector] = contribution
    end
    return expressions
end

function _assemble_frequency_expression!(
    out::CanonicalFrequencyCollision{D,S,O,E1,E2,G,Ctx},
    expression::SpectralDispersiveExpression{C,S,E1,E2,G,Ctx},
    include_external::Bool,
) where {C<:Number,D<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    target = target_family(out)
    parameter = parameters(out)
    for (term, source_coefficient) in expression
        statistical = _collision_statistical_monomial(term, target, include_external)
        source = convert(D, source_coefficient)
        if _has_trotter_shift(term)
            push!(
                out.shifted,
                ShiftedFrequencyCollisionTerm{D,S,E1,E2}(
                    term, source, statistical, parameter
                ),
            )
            continue
        end

        sector = _canonical_frequency_sector(term, statistical, parameter)
        causal = _causal_frequency_expression(term, source, target)
        _merge_causal_frequency_expression!(out.expressions, sector, causal)
    end
    return out
end

"""
    canonical_frequency_collision(collision)

Assemble the complete Kadanoff--Baym collision into canonical affine causal frequency
expressions before any contour integration.

The exact internal statistical factors are part of the sector identity. Terms from the
Kadanoff--Baym distribution coefficient acquire exactly one external `F_target(k)` factor before
aggregation. Consequently cancellations involving routed statistical identity are visible to
the frequency reducer rather than postponed until after termwise integration.

Unshifted spectral/dispersive factors are expanded only by the exact identities
`D=(Gᴿ+Gᴬ)/2` and `A=i(Gᴿ-Gᴬ)`. Finite-Trotter shifted terms remain explicit for the later
regulator layer. No shell delta, principal value, occupation substitution, loop-momentum
quotient, or topology-specific rule is introduced here.
"""
function canonical_frequency_collision(
    collision::SpectralDispersiveCollision{C,S,O,E1,E2,G,Ctx}
) where {C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    D = _canonical_frequency_coefficient_type(C)
    out = CanonicalFrequencyCollision{D,S,O,E1,E2,G,Ctx}(
        Dict{CanonicalFrequencySector{S},CausalFrequencyExpression{D,S}}(),
        ShiftedFrequencyCollisionTerm{D,S,E1,E2}[],
        target_family(collision),
        parameters(collision),
        wigner_context(collision),
    )
    _assemble_frequency_expression!(out, collision_offset(collision), false)
    _assemble_frequency_expression!(
        out, collision_distribution_coefficient(collision), true
    )
    return out
end
