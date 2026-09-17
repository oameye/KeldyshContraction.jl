"""
Typed occupation background used to evaluate a collision linearization.

The evaluator receives an `OccupationAtom` and must return a value convertible to `V`. Declaring
`V` explicitly keeps background evaluation type-stable while allowing numerical or symbolic
number types.
"""
struct OccupationBackground{V<:Number,F}
    evaluator::F
end

function OccupationBackground(::Type{V}, evaluator::F) where {V<:Number,F}
    return OccupationBackground{V,F}(evaluator)
end

function OccupationBackground(evaluator::F, ::Type{V}) where {V<:Number,F}
    return OccupationBackground{V,F}(evaluator)
end

"""Evaluate one occupation atom in a typed background state."""
function occupation_background_value(
    background::OccupationBackground{V}, atom::OccupationAtom
)::V where {V<:Number}
    return convert(V, background.evaluator(atom))
end

"""Evaluate an exact occupation polynomial in a supplied typed background state."""
function evaluate_occupation_polynomial(
    polynomial::OccupationPolynomial{C,S}, background::OccupationBackground{V}
) where {C<:Number,S<:Statistics,V<:Number}
    D = promote_type(C, V)
    value = zero(D)
    for (monomial, coefficient) in polynomial
        term = convert(D, coefficient)
        for atom in monomial
            term *= convert(D, occupation_background_value(background, atom))
        end
        value += term
    end
    return value
end

"""
Occupation-space response coefficients evaluated at a supplied background state.

Each entry retains the exact variation channel `δn_a` while replacing its background occupation
polynomial by the corresponding value at `n̄`.
"""
struct BackgroundOccupationLinearization{C<:Number,S<:Statistics}
    terms::Vector{Pair{OccupationAtom{S},C}}

    function BackgroundOccupationLinearization{C,S}(
        terms::Vector{Pair{OccupationAtom{S},C}}, ::Val{:raw}
    ) where {C<:Number,S<:Statistics}
        return new{C,S}(terms)
    end
end

function _canonical_background_occupation_linearization(
    terms::Vector{Pair{OccupationAtom{S},C}}
) where {C<:Number,S<:Statistics}
    isempty(terms) && return BackgroundOccupationLinearization{C,S}(terms, Val(:raw))
    sort!(terms; by=first)
    out = Pair{OccupationAtom{S},C}[]
    sizehint!(out, length(terms))
    for (atom, coefficient) in terms
        iszero(coefficient) && continue
        if !isempty(out) && isequal(first(out[end]), atom)
            combined = last(out[end]) + coefficient
            pop!(out)
            iszero(combined) || push!(out, atom => combined)
        else
            push!(out, atom => coefficient)
        end
    end
    return BackgroundOccupationLinearization{C,S}(out, Val(:raw))
end

function BackgroundOccupationLinearization{C,S}(
    terms::Vector{Pair{OccupationAtom{S},C}}
) where {C<:Number,S<:Statistics}
    return _canonical_background_occupation_linearization(copy(terms))
end

Base.length(linearization::BackgroundOccupationLinearization) = length(linearization.terms)
Base.isempty(linearization::BackgroundOccupationLinearization) = isempty(linearization.terms)
Base.iszero(linearization::BackgroundOccupationLinearization) = isempty(linearization.terms)
Base.iterate(linearization::BackgroundOccupationLinearization) = iterate(linearization.terms)
function Base.iterate(linearization::BackgroundOccupationLinearization, state)
    return iterate(linearization.terms, state)
end
function Base.eltype(::Type{BackgroundOccupationLinearization{C,S}}) where {C,S}
    return Pair{OccupationAtom{S},C}
end
Base.IteratorSize(::Type{<:BackgroundOccupationLinearization}) = Base.HasLength()

"""Return the exact variation-channel coefficients evaluated at the background state."""
background_occupation_linearization_terms(linearization::BackgroundOccupationLinearization) =
    linearization.terms

"""Evaluate an exact occupation-space Fréchet derivative at a supplied background state."""
function evaluate_occupation_linearization(
    linearization::OccupationLinearization{C,S}, background::OccupationBackground{V}
) where {C<:Number,S<:Statistics,V<:Number}
    D = promote_type(C, V)
    terms = Pair{OccupationAtom{S},D}[]
    sizehint!(terms, length(linearization))
    for (atom, polynomial) in linearization
        coefficient = convert(D, evaluate_occupation_polynomial(polynomial, background))
        push!(terms, atom => coefficient)
    end
    return BackgroundOccupationLinearization{D,S}(terms)
end

"""
Regular strict-QP collision response evaluated around a supplied occupation background.

Collision-sector provenance is unchanged. Each sector now carries numerical or symbolic
coefficients multiplying the explicit variation channels `δn_a`.
"""
struct BackgroundLinearizedCollisionKernel{
    C<:Number,S<:Statistics,O,G,Ctx<:AbstractWignerContext
}
    terms::Dict{CollisionKernelSector{S},BackgroundOccupationLinearization{C,S}}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::BackgroundLinearizedCollisionKernel{C,S,O}) where {C,S,O} = O
statistics(::BackgroundLinearizedCollisionKernel{C,S}) where {C,S} = S
target_family(kernel::BackgroundLinearizedCollisionKernel) = kernel.target
parameters(kernel::BackgroundLinearizedCollisionKernel) = kernel.parameter
gradient_order(::BackgroundLinearizedCollisionKernel{C,S,O,G}) where {C,S,O,G} = Val(G)
wigner_context(kernel::BackgroundLinearizedCollisionKernel) = kernel.context

"""Return the canonical sector-to-background-response mapping."""
background_linearized_collision_terms(kernel::BackgroundLinearizedCollisionKernel) = kernel.terms

Base.length(kernel::BackgroundLinearizedCollisionKernel) = length(kernel.terms)
Base.isempty(kernel::BackgroundLinearizedCollisionKernel) = isempty(kernel.terms)

"""Evaluate a symbolic linearized collision kernel at a supplied occupation background."""
function evaluate_collision_background(
    kernel::LinearizedCollisionKernel{C,S,O,G,Ctx}, background::OccupationBackground{V}
) where {C<:Number,S<:Statistics,O,G,Ctx<:AbstractWignerContext,V<:Number}
    D = promote_type(C, V)
    out = Dict{CollisionKernelSector{S},BackgroundOccupationLinearization{D,S}}()
    for (sector, linearization) in linearized_collision_terms(kernel)
        evaluated = evaluate_occupation_linearization(linearization, background)
        iszero(evaluated) || (out[sector] = evaluated)
    end
    return BackgroundLinearizedCollisionKernel{D,S,O,G,Ctx}(
        out, target_family(kernel), parameters(kernel), wigner_context(kernel)
    )
end

"""Linearize a collision kernel and evaluate its coefficients at a supplied background state."""
function linearize_collision_kernel(kernel::CollisionKernel, background::OccupationBackground)
    return evaluate_collision_background(linearize_collision_kernel(kernel), background)
end

"""Attach an explicit moment descriptor to a background-evaluated collision response."""
function project_collision_moment(
    kernel::K, moment::M
) where {K<:BackgroundLinearizedCollisionKernel,M}
    return CollisionMomentProjection{M,K}(moment, kernel)
end

"""Construct the particle-number moment of a background-evaluated collision response."""
function number_moment_projection(kernel::BackgroundLinearizedCollisionKernel)
    return project_collision_moment(kernel, NumberMoment())
end

"""Construct the energy moment of a background-evaluated collision response."""
function energy_moment_projection(kernel::BackgroundLinearizedCollisionKernel)
    return project_collision_moment(kernel, EnergyMoment())
end

"""Linearize and background-evaluate the collision carried by a formal moment projection."""
function linearize_collision_moment(
    projection::CollisionMomentProjection{M,K}, background::OccupationBackground
) where {M,K<:CollisionKernel}
    linearized = linearize_collision_kernel(projected_collision(projection), background)
    return CollisionMomentProjection{M,typeof(linearized)}(
        moment_test_function(projection), linearized
    )
end

"""Background-evaluate an already symbolic-linearized collision moment projection."""
function linearize_collision_moment(
    projection::CollisionMomentProjection{M,K}, background::OccupationBackground
) where {M,K<:LinearizedCollisionKernel}
    linearized = evaluate_collision_background(projected_collision(projection), background)
    return CollisionMomentProjection{M,typeof(linearized)}(
        moment_test_function(projection), linearized
    )
end
