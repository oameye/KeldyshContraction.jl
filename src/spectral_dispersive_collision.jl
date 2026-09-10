"""Exact causal decomposition kind carried by one internal collision factor."""
@enum SpectralDispersiveKind::UInt8 begin
    CollisionDispersive = 0
    CollisionSpectral = 1
end

"""
One graph contribution after exact retarded/advanced decomposition.

`carrier` reuses the exact `KineticTerm` routing, topology, field-family, regularisation,
and momentum-polynomial storage. Its line kinds are normalized to `KineticSpectral` so the
original R/A spelling cannot obstruct duplicate merging. `kinds` records the physical factor
basis: `CollisionDispersive` for `D=(Gᴿ+Gᴬ)/2` and `CollisionSpectral` for
`A=im*(Gᴿ-Gᴬ)`.
"""
struct SpectralDispersiveTerm{S<:Statistics,E1,E2}
    carrier::KineticTerm{S,E1,E2}
    kinds::FixedVector{E1,SpectralDispersiveKind}
end

statistics(::SpectralDispersiveTerm{S}) where {S} = S
"""Return the canonical spectral/dispersive kind vector for one collision term."""
spectral_dispersive_kinds(term::SpectralDispersiveTerm) = term.kinds
kinematic_factor(term::SpectralDispersiveTerm) = kinematic_factor(term.carrier)
momentum_basis(term::SpectralDispersiveTerm) = momentum_basis(term.carrier)
function external_wigner_momentum(term::SpectralDispersiveTerm)
    return external_wigner_momentum(term.carrier)
end
topology(term::SpectralDispersiveTerm) = term.carrier.topology

function Base.isequal(
    a::SpectralDispersiveTerm{S,E1,E2}, b::SpectralDispersiveTerm{S,E1,E2}
) where {S,E1,E2}
    return isequal(a.carrier, b.carrier) && isequal(a.kinds, b.kinds)
end
Base.:(==)(a::SpectralDispersiveTerm, b::SpectralDispersiveTerm) = isequal(a, b)
function Base.hash(term::SpectralDispersiveTerm, h::UInt)
    return hash(term.kinds, hash(term.carrier, hash(SpectralDispersiveTerm, h)))
end

"""Concrete exact polynomial in spectral/dispersive graph terms."""
struct SpectralDispersiveExpression{
    C<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext
}
    terms::Dict{SpectralDispersiveTerm{S,E1,E2},C}
    context::Ctx
end

function SpectralDispersiveExpression{C,S,E1,E2,G,Ctx}(
    context::Ctx
) where {C<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext}
    K = SpectralDispersiveTerm{S,E1,E2}
    return SpectralDispersiveExpression{C,S,E1,E2,G,Ctx}(Dict{K,C}(), context)
end

statistics(::SpectralDispersiveExpression{C,S}) where {C,S} = S
gradient_order(::SpectralDispersiveExpression{C,S,E1,E2,G}) where {C,S,E1,E2,G} = Val(G)
wigner_context(expression::SpectralDispersiveExpression) = expression.context
Base.length(expression::SpectralDispersiveExpression) = length(expression.terms)
Base.isempty(expression::SpectralDispersiveExpression) = isempty(expression.terms)
Base.iszero(expression::SpectralDispersiveExpression) = isempty(expression.terms)
Base.iterate(expression::SpectralDispersiveExpression) = iterate(expression.terms)
function Base.iterate(expression::SpectralDispersiveExpression, state)
    return iterate(expression.terms, state)
end
function Base.eltype(
    ::Type{SpectralDispersiveExpression{C,S,E1,E2,G,Ctx}}
) where {C,S,E1,E2,G,Ctx}
    return Pair{SpectralDispersiveTerm{S,E1,E2},C}
end

function Base.isequal(
    a::SpectralDispersiveExpression{C,S,E1,E2,G,Ctx},
    b::SpectralDispersiveExpression{C,S,E1,E2,G,Ctx},
) where {C,S,E1,E2,G,Ctx}
    return isequal(a.terms, b.terms) && isequal(a.context, b.context)
end
Base.:(==)(a::SpectralDispersiveExpression, b::SpectralDispersiveExpression) = isequal(a, b)
function Base.hash(expression::SpectralDispersiveExpression, h::UInt)
    return hash(
        expression.context, hash(expression.terms, hash(SpectralDispersiveExpression, h))
    )
end

function Base.push!(
    expression::SpectralDispersiveExpression{C,S,E1,E2,G,Ctx},
    term::SpectralDispersiveTerm{S,E1,E2},
    coefficient::Number,
) where {C<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext}
    value = _simplify(convert(C, coefficient))
    iszero(value) && return expression
    if haskey(expression.terms, term)
        combined = _simplify(expression.terms[term] + value)
        if iszero(combined)
            delete!(expression.terms, term)
        else
            expression.terms[term] = combined
        end
    else
        expression.terms[term] = value
    end
    return expression
end

"""
Exact off-shell collision expression in the dispersive/spectral line basis.

No shell delta function, principal-value denominator, or occupation approximation is encoded
by this type.
"""
struct SpectralDispersiveCollision{
    C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext
}
    offset::SpectralDispersiveExpression{C,S,E1,E2,G,Ctx}
    distribution_coefficient::SpectralDispersiveExpression{C,S,E1,E2,G,Ctx}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::SpectralDispersiveCollision{C,S,O}) where {C,S,O} = O
statistics(::SpectralDispersiveCollision{C,S}) where {C,S} = S
parameters(collision::SpectralDispersiveCollision) = collision.parameter
gradient_order(::SpectralDispersiveCollision{C,S,O,E1,E2,G}) where {C,S,O,E1,E2,G} = Val(G)
wigner_context(collision::SpectralDispersiveCollision) = collision.context
target_family(collision::SpectralDispersiveCollision) = collision.target
collision_offset(collision::SpectralDispersiveCollision) = collision.offset
function collision_distribution_coefficient(collision::SpectralDispersiveCollision)
    return collision.distribution_coefficient
end

function Base.isequal(
    a::SpectralDispersiveCollision{C,S,O,E1,E2,G,Ctx},
    b::SpectralDispersiveCollision{C,S,O,E1,E2,G,Ctx},
) where {C,S,O,E1,E2,G,Ctx}
    return isequal(a.offset, b.offset) &&
           isequal(a.distribution_coefficient, b.distribution_coefficient) &&
           isequal(a.target, b.target) &&
           isequal(a.parameter, b.parameter) &&
           isequal(a.context, b.context)
end
Base.:(==)(a::SpectralDispersiveCollision, b::SpectralDispersiveCollision) = isequal(a, b)
function Base.hash(collision::SpectralDispersiveCollision, h::UInt)
    return hash(
        (
            collision.offset,
            collision.distribution_coefficient,
            collision.target,
            collision.parameter,
            collision.context,
        ),
        hash(SpectralDispersiveCollision, h),
    )
end

@inline function _spectral_dispersive_components(
    line::KineticLine{S}, ::Type{C}
) where {C<:Number,S<:Statistics}
    kind = kinetic_line_kind(line)
    if kind === KineticSpectral
        return Pair{SpectralDispersiveKind,C}[CollisionSpectral => one(C)]
    end

    statistical_weight(line) === NoStatisticalWeight ||
        error("causal kinetic line unexpectedly carries a statistical distribution weight")
    half_i = convert(C, (1 // 2) * im)
    spectral_coefficient = if kind === KineticRetarded
        -half_i
    elseif kind === KineticAdvanced
        half_i
    else
        error("unsupported kinetic line kind in causal decomposition")
    end
    return Pair{SpectralDispersiveKind,C}[
        CollisionDispersive => one(C), CollisionSpectral => spectral_coefficient
    ]
end

@inline function _normalized_collision_line(line::KineticLine{S}) where {S<:Statistics}
    return KineticLine{S}(
        line.family,
        line.out_position,
        line.in_position,
        line.regularisation_shift,
        momentum(line),
        KineticSpectral,
        statistical_weight(line),
    )
end

@inline function _spectral_dispersive_pair_isless(a, b)
    line_a, kind_a = a
    line_b, kind_b = b
    isequal(line_a, line_b) || return isless(line_a, line_b)
    return Int(kind_a) < Int(kind_b)
end

function _spectral_dispersive_term(
    source::KineticTerm{S,E1,E2}, kinds::Vector{SpectralDispersiveKind}
) where {S<:Statistics,E1,E2}
    length(kinds) == E1 || error("causal decomposition changed kinetic line count")
    pairs = Tuple{KineticLine{S},SpectralDispersiveKind}[
        (_normalized_collision_line(line), kind) for
        (line, kind) in zip(kinetic_lines(source), kinds)
    ]
    sort!(pairs; lt=_spectral_dispersive_pair_isless)
    lines = KineticLine{S}[first(pair) for pair in pairs]
    canonical_kinds = SpectralDispersiveKind[last(pair) for pair in pairs]
    carrier = KineticTerm{S,E1,E2}(
        KineticMonomial(lines, Val(E1)),
        source.topology,
        source.basis,
        source.external_momentum,
        source.kinematic,
    )
    return SpectralDispersiveTerm{S,E1,E2}(
        carrier, FixedVector{E1,SpectralDispersiveKind}(canonical_kinds)
    )
end

function _spectral_dispersive_expression(
    expression::KineticExpression{C,S,E1,E2,G,Ctx}
) where {C<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext}
    out = SpectralDispersiveExpression{C,S,E1,E2,G,Ctx}(wigner_context(expression))
    for (source_term, source_coefficient) in expression
        partials = Pair{Vector{SpectralDispersiveKind},C}[SpectralDispersiveKind[] => source_coefficient]
        for line in kinetic_lines(source_term)
            components = _spectral_dispersive_components(line, C)
            next = Pair{Vector{SpectralDispersiveKind},C}[]
            sizehint!(next, length(partials) * length(components))
            for (partial_kinds, partial_coefficient) in partials
                for (kind, component_coefficient) in components
                    coefficient = partial_coefficient * component_coefficient
                    iszero(coefficient) && continue
                    kinds = copy(partial_kinds)
                    push!(kinds, kind)
                    push!(next, kinds => coefficient)
                end
            end
            partials = next
        end
        for (kinds, coefficient) in partials
            push!(out, _spectral_dispersive_term(source_term, kinds), coefficient)
        end
    end
    return out
end

"""
    spectral_dispersive_collision(collision::OffShellCollisionExpression)

Apply the exact identities `Gᴿ=D-im*A/2` and `Gᴬ=D+im*A/2` to the complete off-shell
collision expression. Keldysh-derived spectral factors remain spectral and retain their
statistical weight.

This transformation performs no quasiparticle, shell, delta-function, principal-value, or
occupation approximation. Exact routing, equal-time regularisation shifts, derivative
kinematics, parameter provenance, and the external target family are preserved.
"""
function spectral_dispersive_collision(
    collision::OffShellCollisionExpression{C,S,O,E1,E2,G,Ctx}
) where {C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    return SpectralDispersiveCollision{C,S,O,E1,E2,G,Ctx}(
        _spectral_dispersive_expression(collision.offset),
        _spectral_dispersive_expression(collision.distribution_coefficient),
        collision.target,
        collision.parameter,
        collision.context,
    )
end
