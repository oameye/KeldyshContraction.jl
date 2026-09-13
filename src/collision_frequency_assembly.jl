"""
Physical regular collision sector after exact loop-frequency reduction.

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
Canonical collision-level identity of one linearly dependent shell sector.

The singular affine shell geometry and its finite residual support are kept separate. As for the
regular sector, coordinate topology is deliberately absent so physically identical singular
terms can cancel after their statistical factors have been canonicalized by routed momentum.
"""
struct ReducedDependentCollisionSector{S<:Statistics}
    parameter::ParameterMonomial
    basis::MomentumBasis
    external_momentum::MomentumVariable
    kinematic::MomentumPolynomial{ComplexRationals}
    dependent_support::DependentShellSupport{S}
    residual_support::FrequencySupport{S}
end

statistics(::ReducedDependentCollisionSector{S}) where {S<:Statistics} = S
parameters(sector::ReducedDependentCollisionSector) = sector.parameter
momentum_basis(sector::ReducedDependentCollisionSector) = sector.basis
external_wigner_momentum(sector::ReducedDependentCollisionSector) = sector.external_momentum
kinematic_factor(sector::ReducedDependentCollisionSector) = sector.kinematic
dependent_shell_support(sector::ReducedDependentCollisionSector) = sector.dependent_support
frequency_support(sector::ReducedDependentCollisionSector) = sector.residual_support

function Base.isequal(
    a::ReducedDependentCollisionSector{S}, b::ReducedDependentCollisionSector{S}
) where {S}
    return isequal(a.parameter, b.parameter) &&
           isequal(a.basis, b.basis) &&
           isequal(a.external_momentum, b.external_momentum) &&
           isequal(a.kinematic, b.kinematic) &&
           isequal(a.dependent_support, b.dependent_support) &&
           isequal(a.residual_support, b.residual_support)
end
function Base.:(==)(a::ReducedDependentCollisionSector, b::ReducedDependentCollisionSector)
    return isequal(a, b)
end
function Base.hash(sector::ReducedDependentCollisionSector, h::UInt)
    h = hash(ReducedDependentCollisionSector, h)
    h = hash(sector.parameter, h)
    h = hash(sector.basis, h)
    h = hash(sector.external_momentum, h)
    h = hash(sector.kinematic, h)
    h = hash(sector.dependent_support, h)
    return hash(sector.residual_support, h)
end

"""
Canonical collision-level identity of one unresolved causal frequency sector.

The numerical coefficient of the causal product is deliberately not part of this key: it is
multiplied into the canonical statistical polynomial so source terms with identical physical
singular geometry can merge or cancel exactly. `active_loop_basis_indices` identifies the
original momentum-basis variables represented by the causal denominator coordinates. Topology
and source line indices are provenance only and therefore absent from the physical identity.
"""
struct ReducedCausalCollisionSector{S<:Statistics}
    parameter::ParameterMonomial
    basis::MomentumBasis
    external_momentum::MomentumVariable
    kinematic::MomentumPolynomial{ComplexRationals}
    dependent_support::DependentShellSupport{S}
    support::FrequencySupport{S}
    active_loop_basis_indices::Vector{Int}
    kind::CausalExceptionalKind
    denominators::Vector{CausalFrequencyDenominator{S}}
end

statistics(::ReducedCausalCollisionSector{S}) where {S<:Statistics} = S
parameters(sector::ReducedCausalCollisionSector) = sector.parameter
momentum_basis(sector::ReducedCausalCollisionSector) = sector.basis
external_wigner_momentum(sector::ReducedCausalCollisionSector) = sector.external_momentum
kinematic_factor(sector::ReducedCausalCollisionSector) = sector.kinematic
dependent_shell_support(sector::ReducedCausalCollisionSector) = sector.dependent_support
frequency_support(sector::ReducedCausalCollisionSector) = sector.support
function active_loop_basis_indices(sector::ReducedCausalCollisionSector)
    return sector.active_loop_basis_indices
end
causal_exceptional_kind(sector::ReducedCausalCollisionSector) = sector.kind
causal_denominators(sector::ReducedCausalCollisionSector) = sector.denominators

function Base.isequal(
    a::ReducedCausalCollisionSector{S}, b::ReducedCausalCollisionSector{S}
) where {S}
    return isequal(a.parameter, b.parameter) &&
           isequal(a.basis, b.basis) &&
           isequal(a.external_momentum, b.external_momentum) &&
           isequal(a.kinematic, b.kinematic) &&
           isequal(a.dependent_support, b.dependent_support) &&
           isequal(a.support, b.support) &&
           isequal(a.active_loop_basis_indices, b.active_loop_basis_indices) &&
           a.kind === b.kind &&
           isequal(a.denominators, b.denominators)
end
Base.:(==)(a::ReducedCausalCollisionSector, b::ReducedCausalCollisionSector) = isequal(a, b)
function Base.hash(sector::ReducedCausalCollisionSector, h::UInt)
    h = hash(ReducedCausalCollisionSector, h)
    h = hash(sector.parameter, h)
    h = hash(sector.basis, h)
    h = hash(sector.external_momentum, h)
    h = hash(sector.kinematic, h)
    h = hash(sector.dependent_support, h)
    h = hash(sector.support, h)
    h = hash(sector.active_loop_basis_indices, h)
    h = hash(sector.kind, h)
    return hash(sector.denominators, h)
end

"""Collision provenance multiplying one still-unresolved shifted frequency state."""
struct ReducedTrotterCollisionTerm{C<:Number,S<:Statistics,E1,E2}
    branch::FrequencyIntegrationState{S,E1,E2}
    source_coefficient::C
    statistical::StatisticalMonomial{S}
    parameter::ParameterMonomial
end

source_coefficient(term::ReducedTrotterCollisionTerm) = term.source_coefficient
statistical_monomial(term::ReducedTrotterCollisionTerm) = term.statistical
parameters(term::ReducedTrotterCollisionTerm) = term.parameter
frequency_branch(term::ReducedTrotterCollisionTerm) = term.branch

"""
Complete collision-level result of exact loop-frequency reduction.

Regular, dependent-shell, and unresolved causal sectors are assembled across the offset and
external-distribution pieces of the Kadanoff--Baym identity. Their values are canonical
statistical polynomials, so physically identical source terms can cancel before any occupation
projection. Still-unresolved Trotter states retain their source provenance and are never assigned
a strict-quasiparticle occupation value.
"""
struct ReducedFrequencyCollision{
    C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext
}
    regular::Dict{ReducedCollisionSector{S},StatisticalPolynomial{C,S}}
    dependent::Dict{ReducedDependentCollisionSector{S},StatisticalPolynomial{C,S}}
    causal::Dict{ReducedCausalCollisionSector{S},StatisticalPolynomial{C,S}}
    trotter::Vector{ReducedTrotterCollisionTerm{C,S,E1,E2}}
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
reduced_dependent_terms(collision::ReducedFrequencyCollision) = collision.dependent
reduced_causal_terms(collision::ReducedFrequencyCollision) = collision.causal
reduced_trotter_terms(collision::ReducedFrequencyCollision) = collision.trotter

function Base.isempty(collision::ReducedFrequencyCollision)
    return isempty(collision.regular) &&
           isempty(collision.dependent) &&
           isempty(collision.causal) &&
           isempty(collision.trotter)
end

@inline function _collision_reduction_coefficient_type(::Type{C}) where {C<:Number}
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
    monomial = _internal_statistical_monomial(term)
    include_external || return monomial
    external = StatisticalMonomial(
        StatisticalAtom{S}[_external_statistical_atom(term, target)]
    )
    return monomial * external
end

function _reduced_collision_sector(
    term::SpectralDispersiveTerm{S},
    support::FrequencySupport{S},
    parameter::ParameterMonomial,
) where {S<:Statistics}
    return ReducedCollisionSector{S}(
        parameter,
        momentum_basis(term),
        external_wigner_momentum(term),
        kinematic_factor(term),
        support,
    )
end

function _reduced_dependent_collision_sector(
    term::SpectralDispersiveTerm{S},
    branch::DependentFrequencyTerm{S},
    residual_support::FrequencySupport{S},
    parameter::ParameterMonomial,
) where {S<:Statistics}
    return ReducedDependentCollisionSector{S}(
        parameter,
        momentum_basis(term),
        external_wigner_momentum(term),
        kinematic_factor(term),
        dependent_shell_support(branch),
        residual_support,
    )
end

function _reduced_causal_collision_sector(
    term::SpectralDispersiveTerm{S},
    branch::CausalExceptionalFrequencyTerm{S},
    parameter::ParameterMonomial,
) where {S<:Statistics}
    return ReducedCausalCollisionSector{S}(
        parameter,
        momentum_basis(term),
        external_wigner_momentum(term),
        kinematic_factor(term),
        dependent_shell_support(branch),
        branch.support,
        copy(branch.state.active_loop_basis_indices),
        branch.kind,
        copy(branch.causal_term.denominators),
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

function _reduce_collision_expression!(
    out::ReducedFrequencyCollision{D,S,O,E1,E2,G,Ctx},
    expression::SpectralDispersiveExpression{C,S,E1,E2,G,Ctx},
    include_external::Bool,
) where {C<:Number,D<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    target = target_family(out)
    parameter = parameters(out)
    for (term, source_coefficient) in expression
        statistical = _collision_statistical_monomial(term, target, include_external)
        source = convert(D, source_coefficient)
        result = reduce_frequency_term(term, target)

        for (support, frequency_coefficient) in result.regular
            coefficient = source * convert(D, frequency_coefficient)
            sector = _reduced_collision_sector(term, support, parameter)
            _push_reduced_polynomial!(out.regular, sector, statistical, coefficient)
        end
        for branch in result.dependent
            for (residual_support, frequency_coefficient) in
                dependent_residual_support(branch)
                coefficient = source * convert(D, frequency_coefficient)
                sector = _reduced_dependent_collision_sector(
                    term, branch, residual_support, parameter
                )
                _push_reduced_polynomial!(out.dependent, sector, statistical, coefficient)
            end
        end
        for branch in result.causal
            coefficient = source * convert(D, branch.causal_term.coefficient)
            sector = _reduced_causal_collision_sector(term, branch, parameter)
            _push_reduced_polynomial!(out.causal, sector, statistical, coefficient)
        end
        for branch in result.trotter
            push!(
                out.trotter,
                ReducedTrotterCollisionTerm{D,S,E1,E2}(
                    branch, source, statistical, parameter
                ),
            )
        end
    end
    return out
end

"""
    reduce_frequency_collision(collision)

Apply the exact termwise frequency reducer to the complete Kadanoff--Baym collision expression
and assemble its physical collision-level provenance. Offset terms carry only their internal
statistical factors. Distribution-coefficient terms acquire exactly one external
`F_target(k)` factor before canonical aggregation.

Regular, dependent-shell, and causal singular terms are canonicalized without topology labels.
This permits exact cross-diagram cancellation while preserving genuinely singular support as a
separate category. No occupation substitution, loop-momentum quotient, or finite prescription
for dependent/pinch support is performed here.
"""
function reduce_frequency_collision(
    collision::SpectralDispersiveCollision{C,S,O,E1,E2,G,Ctx}
) where {C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    D = _collision_reduction_coefficient_type(C)
    out = ReducedFrequencyCollision{D,S,O,E1,E2,G,Ctx}(
        Dict{ReducedCollisionSector{S},StatisticalPolynomial{D,S}}(),
        Dict{ReducedDependentCollisionSector{S},StatisticalPolynomial{D,S}}(),
        Dict{ReducedCausalCollisionSector{S},StatisticalPolynomial{D,S}}(),
        ReducedTrotterCollisionTerm{D,S,E1,E2}[],
        target_family(collision),
        parameters(collision),
        wigner_context(collision),
    )
    _reduce_collision_expression!(out, collision_offset(collision), false)
    _reduce_collision_expression!(out, collision_distribution_coefficient(collision), true)
    return out
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

Project only the regular strict-quasiparticle branch to occupation variables using
`F_S = 1 + 2σ_S n_S` and `C_n = σ_S I_package/2`. Singular dependent, causal, and Trotter
branches remain in the input `ReducedFrequencyCollision` and are deliberately absent from this
regular occupation expression.
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
