"""One axis component of an exact routed momentum."""
struct MomentumComponent
    momentum::LinearMomentum
    axis::Symbol

    function MomentumComponent(momentum::LinearMomentum, axis::Symbol)
        derivative_axis_slot(axis)
        return new(momentum, axis)
    end
end

momentum(component::MomentumComponent) = component.momentum
momentum_axis(component::MomentumComponent) = component.axis

function Base.isequal(a::MomentumComponent, b::MomentumComponent)
    return a.axis === b.axis && isequal(a.momentum, b.momentum)
end
Base.:(==)(a::MomentumComponent, b::MomentumComponent) = isequal(a, b)
function Base.hash(component::MomentumComponent, h::UInt)
    return hash(MomentumComponent, hash(component.axis, hash(component.momentum, h)))
end

@inline function _linear_momentum_isless(a::LinearMomentum, b::LinearMomentum)
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        a[i] == b[i] || return a[i] < b[i]
    end
    return length(a) < length(b)
end

function Base.isless(a::MomentumComponent, b::MomentumComponent)
    a_slot = derivative_axis_slot(a.axis)
    b_slot = derivative_axis_slot(b.axis)
    a_slot == b_slot || return a_slot < b_slot
    isequal(a.momentum, b.momentum) && return false
    return _linear_momentum_isless(a.momentum, b.momentum)
end

"""Canonical commutative product of routed momentum components."""
struct MomentumMonomial
    factors::Vector{MomentumComponent}

    function MomentumMonomial(factors::Vector{MomentumComponent}, ::Val{:raw})
        return new(factors)
    end
end

MomentumMonomial() = MomentumMonomial(MomentumComponent[], Val(:raw))
function MomentumMonomial(factors::AbstractVector{MomentumComponent})
    canonical = collect(factors)
    sort!(canonical)
    return MomentumMonomial(canonical, Val(:raw))
end

Base.length(monomial::MomentumMonomial) = length(monomial.factors)
Base.isempty(monomial::MomentumMonomial) = isempty(monomial.factors)
Base.iterate(monomial::MomentumMonomial) = iterate(monomial.factors)
Base.iterate(monomial::MomentumMonomial, state) = iterate(monomial.factors, state)
Base.eltype(::Type{MomentumMonomial}) = MomentumComponent
Base.IteratorSize(::Type{MomentumMonomial}) = Base.HasLength()
Base.isequal(a::MomentumMonomial, b::MomentumMonomial) = isequal(a.factors, b.factors)
Base.:(==)(a::MomentumMonomial, b::MomentumMonomial) = isequal(a, b)
function Base.hash(monomial::MomentumMonomial, h::UInt)
    return hash(MomentumMonomial, hash(monomial.factors, h))
end

function Base.isless(a::MomentumMonomial, b::MomentumMonomial)
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        ai = a.factors[i]
        bi = b.factors[i]
        isequal(ai, bi) || return isless(ai, bi)
    end
    return length(a) < length(b)
end

function Base.:*(a::MomentumMonomial, b::MomentumMonomial)
    return MomentumMonomial(vcat(a.factors, b.factors))
end

"""
Canonical finite polynomial in exact routed momentum components.

The coefficient belongs to the kinematic polynomial, not to the numerical diagram
coefficient. Fourier derivative lowering uses `ComplexRationals` coefficients so all
powers of `i` and integer/rational prefactors remain exact.
"""
struct MomentumPolynomial{C<:Number}
    terms::Vector{Pair{MomentumMonomial,C}}

    function MomentumPolynomial{C}(
        terms::Vector{Pair{MomentumMonomial,C}}, ::Val{:raw}
    ) where {C<:Number}
        return new{C}(terms)
    end
end

function MomentumPolynomial{C}() where {C<:Number}
    return MomentumPolynomial{C}(Pair{MomentumMonomial,C}[], Val(:raw))
end

function _canonical_momentum_polynomial(
    terms::Vector{Pair{MomentumMonomial,C}}
) where {C<:Number}
    isempty(terms) && return MomentumPolynomial{C}()
    sort!(terms; by=first)

    out = Pair{MomentumMonomial,C}[]
    sizehint!(out, length(terms))
    for (monomial, coefficient) in terms
        iszero(coefficient) && continue
        if !isempty(out) && isequal(first(out[end]), monomial)
            combined = last(out[end]) + coefficient
            pop!(out)
            iszero(combined) || push!(out, monomial => combined)
        else
            push!(out, monomial => coefficient)
        end
    end
    return MomentumPolynomial{C}(out, Val(:raw))
end

function MomentumPolynomial{C}(terms::Vector{Pair{MomentumMonomial,C}}) where {C<:Number}
    return _canonical_momentum_polynomial(copy(terms))
end

function MomentumPolynomial(monomial::MomentumMonomial, coefficient::C) where {C<:Number}
    return if iszero(coefficient)
        MomentumPolynomial{C}()
    else
        MomentumPolynomial{C}([monomial => coefficient])
    end
end

Base.length(polynomial::MomentumPolynomial) = length(polynomial.terms)
Base.isempty(polynomial::MomentumPolynomial) = isempty(polynomial.terms)
Base.iszero(polynomial::MomentumPolynomial) = isempty(polynomial.terms)
Base.iterate(polynomial::MomentumPolynomial) = iterate(polynomial.terms)
Base.iterate(polynomial::MomentumPolynomial, state) = iterate(polynomial.terms, state)
Base.eltype(::Type{MomentumPolynomial{C}}) where {C<:Number} = Pair{MomentumMonomial,C}
Base.IteratorSize(::Type{<:MomentumPolynomial}) = Base.HasLength()
function Base.isequal(a::MomentumPolynomial, b::MomentumPolynomial)
    return isequal(a.terms, b.terms)
end
Base.:(==)(a::MomentumPolynomial, b::MomentumPolynomial) = isequal(a, b)
function Base.hash(polynomial::MomentumPolynomial, h::UInt)
    return hash(MomentumPolynomial, hash(polynomial.terms, h))
end

function Base.:+(
    a::MomentumPolynomial{C1}, b::MomentumPolynomial{C2}
) where {C1<:Number,C2<:Number}
    D = promote_type(C1, C2)
    terms = Pair{MomentumMonomial,D}[]
    sizehint!(terms, length(a) + length(b))
    for (monomial, coefficient) in a
        push!(terms, monomial => convert(D, coefficient))
    end
    for (monomial, coefficient) in b
        push!(terms, monomial => convert(D, coefficient))
    end
    return _canonical_momentum_polynomial(terms)
end

function Base.:-(polynomial::MomentumPolynomial{C}) where {C<:Number}
    terms = Pair{MomentumMonomial,C}[
        monomial => -coefficient for (monomial, coefficient) in polynomial
    ]
    return MomentumPolynomial{C}(terms)
end
Base.:-(a::MomentumPolynomial, b::MomentumPolynomial) = a + (-b)

function Base.:*(
    a::MomentumPolynomial{C1}, b::MomentumPolynomial{C2}
) where {C1<:Number,C2<:Number}
    D = promote_type(C1, C2)
    terms = Pair{MomentumMonomial,D}[]
    sizehint!(terms, length(a) * length(b))
    for (ma, ca) in a, (mb, cb) in b
        push!(terms, ma * mb => convert(D, ca) * convert(D, cb))
    end
    return _canonical_momentum_polynomial(terms)
end

function Base.:*(
    coefficient::D, polynomial::MomentumPolynomial{C}
) where {D<:Number,C<:Number}
    P = promote_type(C, D)
    terms = Pair{MomentumMonomial,P}[
        monomial => convert(P, coefficient) * convert(P, value) for
        (monomial, value) in polynomial
    ]
    return _canonical_momentum_polynomial(terms)
end
Base.:*(polynomial::MomentumPolynomial, coefficient::Number) = coefficient * polynomial

function _kinematic_identity()
    return MomentumPolynomial(MomentumMonomial(), one(ComplexRationals))
end

@inline function _fourier_axis_metric(axis::Symbol)::Int8
    axis === :t && return Int8(-1)
    derivative_axis_slot(axis)
    return Int8(1)
end

@inline function _fourier_derivative_phase(axis::Symbol, outgoing::Bool)::ComplexRationals
    endpoint_sign = outgoing ? Int8(1) : Int8(-1)
    sign = endpoint_sign * _fourier_axis_metric(axis)
    return ComplexRationals(0 // 1, Int(sign) // 1)
end

function _derivative_factor(momentum::LinearMomentum, axis::Symbol, outgoing::Bool)
    component = MomentumComponent(momentum, axis)
    monomial = MomentumMonomial(MomentumComponent[component])
    return MomentumPolynomial(monomial, _fourier_derivative_phase(axis, outgoing))
end

"""
Lower coordinate derivatives on a routed diagram to an exact momentum polynomial.

The convention is the homogeneous inverse of
`A(X,k) = ∫ ds exp[-i(k⋅s - εt)] A(X+s/2,X-s/2)`. Therefore spatial
derivatives on an edge's out endpoint contribute `+i p_α`, spatial derivatives on
the in endpoint contribute `-i p_α`, and time derivatives carry the opposite sign.
The exact routed `LinearMomentum` is retained as one factor rather than expanded into
an external/loop symbolic AST.
"""
function lower_fourier_derivatives(
    diagram::FourierDiagram{S,E1,E2,Nothing}
) where {S<:Statistics,E1,E2}
    polynomial = _kinematic_identity()
    for (edge, routed_momentum) in
        zip(contractions(diagram.coordinate), diagram.edge_momenta)
        out_field, in_field = fields(edge)
        for axis in derivative_multiindex(out_field)
            polynomial *= _derivative_factor(routed_momentum, axis, true)
        end
        for axis in derivative_multiindex(in_field)
            polynomial *= _derivative_factor(routed_momentum, axis, false)
        end
    end

    return FourierDiagram{S,E1,E2,MomentumPolynomial{ComplexRationals}}(
        diagram.coordinate,
        diagram.basis,
        diagram.edge_momenta,
        diagram.external_count,
        diagram.loop_count,
        polynomial,
    )
end
