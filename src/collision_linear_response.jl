"""
Canonical Fréchet derivative of an occupation polynomial.

Each entry maps one varied occupation atom `δn_a` to the exact background polynomial
multiplying that variation. Repeated factors are differentiated with their exact integer
multiplicity, and equal variation atoms are merged canonically.
"""
struct OccupationLinearization{C<:Number,S<:Statistics}
    terms::Vector{Pair{OccupationAtom{S},OccupationPolynomial{C,S}}}

    function OccupationLinearization{C,S}(
        terms::Vector{Pair{OccupationAtom{S},OccupationPolynomial{C,S}}}, ::Val{:raw}
    ) where {C<:Number,S<:Statistics}
        return new{C,S}(terms)
    end
end

function _canonical_occupation_linearization(
    terms::Vector{Pair{OccupationAtom{S},OccupationPolynomial{C,S}}}
) where {C<:Number,S<:Statistics}
    isempty(terms) && return OccupationLinearization{C,S}(terms, Val(:raw))
    sort!(terms; by=first)
    out = Pair{OccupationAtom{S},OccupationPolynomial{C,S}}[]
    sizehint!(out, length(terms))
    for (atom, polynomial) in terms
        iszero(polynomial) && continue
        if !isempty(out) && isequal(first(out[end]), atom)
            combined = last(out[end]) + polynomial
            pop!(out)
            iszero(combined) || push!(out, atom => combined)
        else
            push!(out, atom => polynomial)
        end
    end
    return OccupationLinearization{C,S}(out, Val(:raw))
end

function OccupationLinearization{C,S}(
    terms::Vector{Pair{OccupationAtom{S},OccupationPolynomial{C,S}}}
) where {C<:Number,S<:Statistics}
    return _canonical_occupation_linearization(copy(terms))
end

Base.length(linearization::OccupationLinearization) = length(linearization.terms)
Base.isempty(linearization::OccupationLinearization) = isempty(linearization.terms)
Base.iszero(linearization::OccupationLinearization) = isempty(linearization.terms)
Base.iterate(linearization::OccupationLinearization) = iterate(linearization.terms)
function Base.iterate(linearization::OccupationLinearization, state)
    return iterate(linearization.terms, state)
end
function Base.eltype(::Type{OccupationLinearization{C,S}}) where {C,S}
    return Pair{OccupationAtom{S},OccupationPolynomial{C,S}}
end
Base.IteratorSize(::Type{<:OccupationLinearization}) = Base.HasLength()

"""Return the canonical variation-atom to background-polynomial entries."""
occupation_linearization_terms(linearization::OccupationLinearization) = linearization.terms

"""
    occupation_linearization(polynomial)

Construct the exact first Fréchet derivative of a canonical occupation polynomial.

For example, `c*n(k)^2*n(q)` contributes `2c*n(k)*n(q)` to the `δn(k)` channel and
`c*n(k)^2` to the `δn(q)` channel. No background state is substituted at this stage.
"""
function occupation_linearization(
    polynomial::OccupationPolynomial{C,S}
) where {C<:Number,S<:Statistics}
    D = promote_type(C, Int)
    raw = Pair{OccupationAtom{S},OccupationPolynomial{D,S}}[]

    for (monomial, coefficient) in polynomial
        factors = monomial.factors
        i = firstindex(factors)
        while i <= lastindex(factors)
            atom = factors[i]
            j = i
            while j < lastindex(factors) && isequal(factors[j + 1], atom)
                j += 1
            end
            multiplicity = j - i + 1
            residual_factors = copy(factors)
            deleteat!(residual_factors, i)
            residual = OccupationMonomial(residual_factors)
            derivative_coefficient = convert(D, coefficient) * convert(D, multiplicity)
            contribution = OccupationPolynomial{D,S}([residual => derivative_coefficient])
            push!(raw, atom => contribution)
            i = j + 1
        end
    end

    return _canonical_occupation_linearization(raw)
end

"""
Exact linearization of a regular strict-QP collision kernel in occupation space.

The phase-space sector, shell/PV support, derivative kinematics, perturbative provenance and
Wigner context are unchanged. Only each nonlinear occupation polynomial is replaced by its exact
Fréchet derivative.
"""
struct LinearizedCollisionKernel{C<:Number,S<:Statistics,O,G,Ctx<:AbstractWignerContext}
    terms::Dict{CollisionKernelSector{S},OccupationLinearization{C,S}}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::LinearizedCollisionKernel{C,S,O}) where {C,S,O} = O
statistics(::LinearizedCollisionKernel{C,S}) where {C,S} = S
target_family(kernel::LinearizedCollisionKernel) = kernel.target
parameters(kernel::LinearizedCollisionKernel) = kernel.parameter
gradient_order(::LinearizedCollisionKernel{C,S,O,G}) where {C,S,O,G} = Val(G)
wigner_context(kernel::LinearizedCollisionKernel) = kernel.context

"""Return the canonical sector-to-Fréchet-derivative mapping of a linearized collision kernel."""
linearized_collision_terms(kernel::LinearizedCollisionKernel) = kernel.terms

Base.length(kernel::LinearizedCollisionKernel) = length(kernel.terms)
Base.isempty(kernel::LinearizedCollisionKernel) = isempty(kernel.terms)

"""Construct the exact occupation-space Fréchet derivative of a `CollisionKernel`."""
function linearize_collision_kernel(
    kernel::CollisionKernel{C,S,O,G,Ctx}
) where {C<:Number,S<:Statistics,O,G,Ctx<:AbstractWignerContext}
    D = promote_type(C, Int)
    out = Dict{CollisionKernelSector{S},OccupationLinearization{D,S}}()
    for (sector, polynomial) in collision_kernel_terms(kernel)
        linearization = occupation_linearization(polynomial)
        iszero(linearization) || (out[sector] = linearization)
    end
    return LinearizedCollisionKernel{D,S,O,G,Ctx}(
        out, target_family(kernel), parameters(kernel), wigner_context(kernel)
    )
end

"""Marker for the particle-number moment, whose collision test function is one."""
struct NumberMoment end

"""Marker for the single-particle-energy moment of the external kinetic leg."""
struct EnergyMoment end

"""
Formal moment projection of a collision object.

`moment` is deliberately parametric: package-provided `NumberMoment` and `EnergyMoment` have exact
weights, while downstream code may attach any concrete, inspectable user moment descriptor without
forcing a trap model, phase-space measure, numerical quadrature or closure prescription into the
kinetic compiler.
"""
struct CollisionMomentProjection{M,K}
    moment::M
    collision::K
end

"""Return the explicit moment/test-function descriptor."""
moment_test_function(projection::CollisionMomentProjection) = projection.moment

"""Return the nonlinear or linearized collision object carried by the projection."""
projected_collision(projection::CollisionMomentProjection) = projection.collision

"""Attach an explicit moment/test-function descriptor to a nonlinear collision kernel."""
function project_collision_moment(kernel::K, moment::M) where {K<:CollisionKernel,M}
    return CollisionMomentProjection{M,K}(moment, kernel)
end

"""Attach an explicit moment/test-function descriptor to a linearized collision kernel."""
function project_collision_moment(
    kernel::K, moment::M
) where {K<:LinearizedCollisionKernel,M}
    return CollisionMomentProjection{M,K}(moment, kernel)
end

"""Construct the formal particle-number collision moment `Ṅ = ∫ C`."""
function number_moment_projection(kernel::CollisionKernel)
    return project_collision_moment(kernel, NumberMoment())
end
function number_moment_projection(kernel::LinearizedCollisionKernel)
    return project_collision_moment(kernel, NumberMoment())
end

"""Construct the formal energy collision moment `Ė = ∫ ε_k C`."""
function energy_moment_projection(kernel::CollisionKernel)
    return project_collision_moment(kernel, EnergyMoment())
end
function energy_moment_projection(kernel::LinearizedCollisionKernel)
    return project_collision_moment(kernel, EnergyMoment())
end

"""Exact test-function weight for the particle-number moment."""
function moment_weight(::NumberMoment, ::CollisionKernelSector, ::FieldFamily)
    return one(EnergyCoefficient)
end

"""Exact external quasiparticle-energy test function for one canonical collision sector."""
function moment_weight(
    ::EnergyMoment, sector::CollisionKernelSector{S}, target::FieldFamily{S}
) where {S<:Statistics}
    basis = momentum_basis(sector)
    external = external_wigner_momentum(sector)
    external_index = _external_basis_index(basis, external)
    return EnergyForm(DispersionAtom(target, basis_momentum(basis, external_index)))
end

"""Resolve the exact particle-number weight for one sector of a projection."""
function moment_weight(
    projection::CollisionMomentProjection{NumberMoment,K}, sector::CollisionKernelSector
) where {K}
    return moment_weight(
        NumberMoment(), sector, target_family(projected_collision(projection))
    )
end

"""Resolve the exact energy weight for one sector of a projection."""
function moment_weight(
    projection::CollisionMomentProjection{EnergyMoment,K}, sector::CollisionKernelSector
) where {K}
    return moment_weight(
        EnergyMoment(), sector, target_family(projected_collision(projection))
    )
end

"""
Linearize the collision operator carried by a formal moment projection.

The moment descriptor is retained exactly, so subsequent phase-space closure may act on the same
test function before or after occupation-space linearization.
"""
function linearize_collision_moment(
    projection::CollisionMomentProjection{M,K}
) where {M,K<:CollisionKernel}
    linearized = linearize_collision_kernel(projected_collision(projection))
    return CollisionMomentProjection{M,typeof(linearized)}(
        moment_test_function(projection), linearized
    )
end

function linearize_collision_moment(
    projection::CollisionMomentProjection{M,K}
) where {M,K<:LinearizedCollisionKernel}
    return projection
end
