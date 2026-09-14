"""
Physical regular collision sector after exact canonical frequency reduction.

Topology is deliberately absent from this identity. Terms with the same perturbation parameter,
momentum basis, external momentum, derivative kinematics, and canonical shell/PV support share one
statistical polynomial and can therefore cancel across diagram labels.
"""
struct ReducedCollisionSector{S<:Statistics}
    parameter::ParameterMonomial
    basis::MomentumBasis
    external_momentum::MomentumVariable
    kinematic::MomentumPolynomial{ComplexRationals}
    support::FrequencySupport{S}
end

statistics(::ReducedCollisionSector{S}) where {S<:Statistics} = S
parameters(sector::ReducedCollisionSector) = sector.parameter
momentum_basis(sector::ReducedCollisionSector) = sector.basis
external_wigner_momentum(sector::ReducedCollisionSector) = sector.external_momentum
kinematic_factor(sector::ReducedCollisionSector) = sector.kinematic
frequency_support(sector::ReducedCollisionSector) = sector.support

function Base.isequal(a::ReducedCollisionSector{S}, b::ReducedCollisionSector{S}) where {S}
    return isequal(a.parameter, b.parameter) &&
           isequal(a.basis, b.basis) &&
           isequal(a.external_momentum, b.external_momentum) &&
           isequal(a.kinematic, b.kinematic) &&
           isequal(a.support, b.support)
end
Base.:(==)(a::ReducedCollisionSector, b::ReducedCollisionSector) = isequal(a, b)
function Base.hash(sector::ReducedCollisionSector, h::UInt)
    h = hash(ReducedCollisionSector, h)
    h = hash(sector.parameter, h)
    h = hash(sector.basis, h)
    h = hash(sector.external_momentum, h)
    h = hash(sector.kinematic, h)
    return hash(sector.support, h)
end

"""
Complete collision-level result after canonical Trotter/frequency reduction and boundary lowering.

`regular` contains only finite shell/PV support that can be projected to occupation variables.
`blocked` keeps genuine causal-frequency blockers unchanged, and `unresolved` keeps shifted states
that did not admit a structural local equal-time reduction. No finite-width or pinch prescription
is introduced here.
"""
struct ReducedFrequencyCollision{
    C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext
}
    regular::Dict{ReducedCollisionSector{S},StatisticalPolynomial{C,S}}
    blocked::Dict{TrotterFrequencyGroupKey{S},TrotterFrequencyBlockedContribution{C,S}}
    unresolved::Vector{TrotterFrequencyState{C,S,E1,E2}}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::ReducedFrequencyCollision{C,S,O}) where {C,S,O} = O
statistics(::ReducedFrequencyCollision{C,S}) where {C,S} = S
target_family(collision::ReducedFrequencyCollision) = collision.target
parameters(collision::ReducedFrequencyCollision) = collision.parameter
gradient_order(::ReducedFrequencyCollision{C,S,O,E1,E2,G}) where {C,S,O,E1,E2,G} = Val(G)
wigner_context(collision::ReducedFrequencyCollision) = collision.context
reduced_regular_terms(collision::ReducedFrequencyCollision) = collision.regular
reduced_blocked_terms(collision::ReducedFrequencyCollision) = collision.blocked
reduced_trotter_terms(collision::ReducedFrequencyCollision) = collision.unresolved

function Base.isempty(collision::ReducedFrequencyCollision)
    return isempty(collision.regular) &&
           isempty(collision.blocked) &&
           isempty(collision.unresolved)
end

function _empty_frequency_support(::Type{S}) where {S<:Statistics}
    return FrequencySupport(EnergyShell{S}[], PrincipalValueSupport{S}[])
end

function _reduced_collision_sector(
    sector::CanonicalFrequencySector{S}, support::FrequencySupport{S}
) where {S<:Statistics}
    return ReducedCollisionSector{S}(
        sector.parameter, sector.basis, sector.external_momentum, sector.kinematic, support
    )
end

function _push_reduced_polynomial!(
    terms::Dict{K,StatisticalPolynomial{C,S}},
    key::K,
    monomial::StatisticalMonomial{S},
    coefficient::C,
) where {K,C<:Number,S<:Statistics}
    iszero(coefficient) && return terms
    contribution = StatisticalPolynomial{C,S}([monomial => coefficient])
    if haskey(terms, key)
        combined = terms[key] + contribution
        if iszero(combined)
            delete!(terms, key)
        else
            terms[key] = combined
        end
    else
        terms[key] = contribution
    end
    return terms
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

function _lower_causal_frequency_term(
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
            ArgumentError(
                "zero causal energy denominator requires singular-support handling"
            ),
        )

        next_states = State[]
        sizehint!(next_states, (iszero(denominator.infinitesimal) ? 1 : 2) * length(states))
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
        _push_frequency_support_coefficient!(
            out, FrequencySupport(shells, principal_values), coefficient
        )
    end
    return out
end

function _lower_causal_frequency_expression(
    expression::CausalFrequencyExpression{C,S}
) where {C<:Number,S<:Statistics}
    D = promote_type(C, ComplexRationals)
    out = Dict{FrequencySupport{S},D}()
    for term in causal_frequency_terms(expression)
        lowered = _lower_causal_frequency_term(term)
        for (support, coefficient) in lowered
            _push_frequency_support_coefficient!(out, support, coefficient)
        end
    end
    return out
end

"""
    reduce_frequency_collision(collision)

Compile one complete canonical collision through the already certified Trotter/frequency reducer.
Scalar constants become regular contributions with empty frequency support. Terminal finite causal
boundary expressions are lowered to exact shell/PV support only after all frequency-stage joins and
canonical cancellations have occurred. Contributions carrying the statistical monomial embedded in
their `CanonicalFrequencySector` are multiplied by the resulting exact frequency coefficient and
assembled into canonical polynomials at collision level. Typed pinch/deferred branches remain
explicit and are never assigned a finite strict-quasiparticle value.
"""
function reduce_frequency_collision(
    collision::CanonicalFrequencyCollision{C,S,O,E1,E2,G,Ctx}
) where {C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    D = promote_type(C, ComplexRationals)
    regular = Dict{ReducedCollisionSector{S},StatisticalPolynomial{D,S}}()
    reduced = reduce_canonical_trotter_frequencies(collision)

    empty_support = _empty_frequency_support(S)
    for (sector, coefficient) in canonical_trotter_constants(reduced)
        key = _reduced_collision_sector(sector, empty_support)
        _push_reduced_polynomial!(
            regular, key, statistical_monomial(sector), convert(D, coefficient)
        )
    end

    for (sector, expression) in canonical_trotter_boundary_expressions(reduced)
        lowered = _lower_causal_frequency_expression(expression)
        for (support, coefficient) in lowered
            key = _reduced_collision_sector(sector, support)
            _push_reduced_polynomial!(
                regular, key, statistical_monomial(sector), convert(D, coefficient)
            )
        end
    end

    return ReducedFrequencyCollision{D,S,O,E1,E2,G,Ctx}(
        regular,
        copy(blocked_trotter_contributions(reduced)),
        copy(unresolved_trotter_states(reduced)),
        target_family(collision),
        parameters(collision),
        wigner_context(collision),
    )
end

function reduce_frequency_collision(collision::SpectralDispersiveCollision)
    return reduce_frequency_collision(canonical_frequency_collision(collision))
end

"""Regular strict-QP collision after exact `F -> n` substitution and collision normalization."""
struct OccupationReducedExpression{C<:Number,S<:Statistics,O,G,Ctx<:AbstractWignerContext}
    terms::Dict{ReducedCollisionSector{S},OccupationPolynomial{C,S}}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::OccupationReducedExpression{C,S,O}) where {C,S,O} = O
statistics(::OccupationReducedExpression{C,S}) where {C,S} = S
target_family(expression::OccupationReducedExpression) = expression.target
parameters(expression::OccupationReducedExpression) = expression.parameter
gradient_order(::OccupationReducedExpression{C,S,O,G}) where {C,S,O,G} = Val(G)
wigner_context(expression::OccupationReducedExpression) = expression.context
occupation_reduced_terms(expression::OccupationReducedExpression) = expression.terms
Base.length(expression::OccupationReducedExpression) = length(expression.terms)
Base.isempty(expression::OccupationReducedExpression) = isempty(expression.terms)

"""
    occupation_reduced_expression(collision)

Project only the finite regular branch to occupation variables using `F_S = 1 + 2σ_S n_S` and
`C_n = σ_S I_package/2`. Typed causal blockers and unresolved Trotter states remain in the input
`ReducedFrequencyCollision` and are deliberately absent from this regular occupation expression.
"""
function occupation_reduced_expression(
    collision::ReducedFrequencyCollision{C,S,O,E1,E2,G,Ctx}
) where {C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    D = promote_type(C, Rational{Int})
    terms = Dict{ReducedCollisionSector{S},OccupationPolynomial{D,S}}()
    for (sector, statistical) in collision.regular
        occupation = occupation_collision_polynomial(statistical)
        iszero(occupation) || (terms[sector] = occupation)
    end
    return OccupationReducedExpression{D,S,O,G,Ctx}(
        terms, target_family(collision), parameters(collision), wigner_context(collision)
    )
end
