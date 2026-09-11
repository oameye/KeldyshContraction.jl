"""One exact statistical distribution factor `F_family(momentum)`."""
struct StatisticalAtom{S<:Statistics}
    family::FieldFamily{S}
    momentum::LinearMomentum
end

statistics(::StatisticalAtom{S}) where {S<:Statistics} = S
statistical_family(atom::StatisticalAtom) = atom.family
momentum(atom::StatisticalAtom) = atom.momentum

function Base.isequal(a::StatisticalAtom{S}, b::StatisticalAtom{S}) where {S<:Statistics}
    return isequal(a.family, b.family) && isequal(a.momentum, b.momentum)
end
function Base.:(==)(a::StatisticalAtom{S}, b::StatisticalAtom{S}) where {S<:Statistics}
    return isequal(a, b)
end
function Base.hash(atom::StatisticalAtom, h::UInt)
    return hash(atom.momentum, hash(atom.family, hash(StatisticalAtom, h)))
end
function Base.isless(a::StatisticalAtom{S}, b::StatisticalAtom{S}) where {S<:Statistics}
    isequal(a.family, b.family) || return isless(a.family, b.family)
    isequal(a.momentum, b.momentum) && return false
    return _linear_momentum_isless(a.momentum, b.momentum)
end

"""Canonical commutative monomial in statistical distribution atoms."""
struct StatisticalMonomial{S<:Statistics}
    factors::Vector{StatisticalAtom{S}}

    function StatisticalMonomial{S}(
        factors::Vector{StatisticalAtom{S}}, ::Val{:raw}
    ) where {S<:Statistics}
        return new{S}(factors)
    end
end

function StatisticalMonomial{S}() where {S<:Statistics}
    return StatisticalMonomial{S}(StatisticalAtom{S}[], Val(:raw))
end

function StatisticalMonomial(
    factors::AbstractVector{StatisticalAtom{S}}
) where {S<:Statistics}
    canonical = collect(factors)
    sort!(canonical)
    return StatisticalMonomial{S}(canonical, Val(:raw))
end

Base.length(monomial::StatisticalMonomial) = length(monomial.factors)
Base.isempty(monomial::StatisticalMonomial) = isempty(monomial.factors)
Base.iterate(monomial::StatisticalMonomial) = iterate(monomial.factors)
Base.iterate(monomial::StatisticalMonomial, state) = iterate(monomial.factors, state)
Base.eltype(::Type{StatisticalMonomial{S}}) where {S} = StatisticalAtom{S}
Base.IteratorSize(::Type{<:StatisticalMonomial}) = Base.HasLength()
function Base.isequal(a::StatisticalMonomial{S}, b::StatisticalMonomial{S}) where {S}
    return isequal(a.factors, b.factors)
end
Base.:(==)(a::StatisticalMonomial{S}, b::StatisticalMonomial{S}) where {S} = isequal(a, b)
function Base.hash(monomial::StatisticalMonomial, h::UInt)
    return hash(monomial.factors, hash(StatisticalMonomial, h))
end
function Base.isless(a::StatisticalMonomial{S}, b::StatisticalMonomial{S}) where {S}
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        ai = a.factors[i]
        bi = b.factors[i]
        isequal(ai, bi) || return isless(ai, bi)
    end
    return length(a) < length(b)
end
function Base.:*(a::StatisticalMonomial{S}, b::StatisticalMonomial{S}) where {S}
    return StatisticalMonomial(vcat(a.factors, b.factors))
end

"""Canonical finite polynomial in exact statistical distribution factors."""
struct StatisticalPolynomial{C<:Number,S<:Statistics}
    terms::Vector{Pair{StatisticalMonomial{S},C}}

    function StatisticalPolynomial{C,S}(
        terms::Vector{Pair{StatisticalMonomial{S},C}}, ::Val{:raw}
    ) where {C<:Number,S<:Statistics}
        return new{C,S}(terms)
    end
end

function StatisticalPolynomial{C,S}() where {C<:Number,S<:Statistics}
    return StatisticalPolynomial{C,S}(Pair{StatisticalMonomial{S},C}[], Val(:raw))
end

function _canonical_statistical_polynomial(
    terms::Vector{Pair{StatisticalMonomial{S},C}}
) where {C<:Number,S<:Statistics}
    isempty(terms) && return StatisticalPolynomial{C,S}()
    sort!(terms; by=first)
    out = Pair{StatisticalMonomial{S},C}[]
    sizehint!(out, length(terms))
    for (monomial, coefficient) in terms
        iszero(coefficient) && continue
        if !isempty(out) && isequal(first(out[end]), monomial)
            combined = _simplify(last(out[end]) + coefficient)
            pop!(out)
            iszero(combined) || push!(out, monomial => combined)
        else
            push!(out, monomial => coefficient)
        end
    end
    return StatisticalPolynomial{C,S}(out, Val(:raw))
end

function StatisticalPolynomial{C,S}(
    terms::Vector{Pair{StatisticalMonomial{S},C}}
) where {C<:Number,S<:Statistics}
    return _canonical_statistical_polynomial(copy(terms))
end

function StatisticalPolynomial(
    atom::StatisticalAtom{S}, coefficient::C
) where {C<:Number,S<:Statistics}
    monomial = StatisticalMonomial(StatisticalAtom{S}[atom])
    return StatisticalPolynomial{C,S}([monomial => coefficient])
end

Base.length(polynomial::StatisticalPolynomial) = length(polynomial.terms)
Base.isempty(polynomial::StatisticalPolynomial) = isempty(polynomial.terms)
Base.iszero(polynomial::StatisticalPolynomial) = isempty(polynomial.terms)
Base.iterate(polynomial::StatisticalPolynomial) = iterate(polynomial.terms)
Base.iterate(polynomial::StatisticalPolynomial, state) = iterate(polynomial.terms, state)
Base.eltype(::Type{StatisticalPolynomial{C,S}}) where {C,S} = Pair{StatisticalMonomial{S},C}
Base.IteratorSize(::Type{<:StatisticalPolynomial}) = Base.HasLength()
function Base.isequal(
    a::StatisticalPolynomial{C,S}, b::StatisticalPolynomial{C,S}
) where {C,S}
    return isequal(a.terms, b.terms)
end
Base.:(==)(a::StatisticalPolynomial, b::StatisticalPolynomial) = isequal(a, b)
function Base.hash(polynomial::StatisticalPolynomial, h::UInt)
    return hash(polynomial.terms, hash(StatisticalPolynomial, h))
end

function Base.zero(::Type{StatisticalPolynomial{C,S}}) where {C<:Number,S<:Statistics}
    return StatisticalPolynomial{C,S}()
end
Base.zero(polynomial::StatisticalPolynomial{C,S}) where {C,S} = zero(typeof(polynomial))
function Base.one(::Type{StatisticalPolynomial{C,S}}) where {C<:Number,S<:Statistics}
    monomial = StatisticalMonomial{S}()
    return StatisticalPolynomial{C,S}([monomial => one(C)])
end
Base.one(polynomial::StatisticalPolynomial{C,S}) where {C,S} = one(typeof(polynomial))

function Base.:+(
    a::StatisticalPolynomial{C1,S}, b::StatisticalPolynomial{C2,S}
) where {C1<:Number,C2<:Number,S<:Statistics}
    D = promote_type(C1, C2)
    terms = Pair{StatisticalMonomial{S},D}[]
    sizehint!(terms, length(a) + length(b))
    for (monomial, coefficient) in a
        push!(terms, monomial => convert(D, coefficient))
    end
    for (monomial, coefficient) in b
        push!(terms, monomial => convert(D, coefficient))
    end
    return _canonical_statistical_polynomial(terms)
end

function Base.:-(polynomial::StatisticalPolynomial{C,S}) where {C<:Number,S<:Statistics}
    return (-one(C)) * polynomial
end
Base.:-(a::StatisticalPolynomial, b::StatisticalPolynomial) = a + (-b)

function Base.:*(
    a::StatisticalPolynomial{C1,S}, b::StatisticalPolynomial{C2,S}
) where {C1<:Number,C2<:Number,S<:Statistics}
    D = promote_type(C1, C2)
    terms = Pair{StatisticalMonomial{S},D}[]
    sizehint!(terms, length(a) * length(b))
    for (ma, ca) in a, (mb, cb) in b
        push!(terms, ma * mb => convert(D, ca) * convert(D, cb))
    end
    return _canonical_statistical_polynomial(terms)
end

function Base.:*(
    coefficient::D, polynomial::StatisticalPolynomial{C,S}
) where {D<:Number,C<:Number,S<:Statistics}
    P = promote_type(C, D)
    terms = Pair{StatisticalMonomial{S},P}[
        monomial => convert(P, coefficient) * convert(P, value) for
        (monomial, value) in polynomial
    ]
    return _canonical_statistical_polynomial(terms)
end
Base.:*(polynomial::StatisticalPolynomial, coefficient::Number) = coefficient * polynomial

"""One exact occupation factor `n_family(momentum)`."""
struct OccupationAtom{S<:Statistics}
    family::FieldFamily{S}
    momentum::LinearMomentum
end

statistics(::OccupationAtom{S}) where {S<:Statistics} = S
occupation_family(atom::OccupationAtom) = atom.family
momentum(atom::OccupationAtom) = atom.momentum

function Base.isequal(a::OccupationAtom{S}, b::OccupationAtom{S}) where {S<:Statistics}
    return isequal(a.family, b.family) && isequal(a.momentum, b.momentum)
end
Base.:(==)(a::OccupationAtom{S}, b::OccupationAtom{S}) where {S<:Statistics} = isequal(a, b)
function Base.hash(atom::OccupationAtom, h::UInt)
    return hash(atom.momentum, hash(atom.family, hash(OccupationAtom, h)))
end
function Base.isless(a::OccupationAtom{S}, b::OccupationAtom{S}) where {S<:Statistics}
    isequal(a.family, b.family) || return isless(a.family, b.family)
    isequal(a.momentum, b.momentum) && return false
    return _linear_momentum_isless(a.momentum, b.momentum)
end

"""Canonical commutative monomial in occupation atoms."""
struct OccupationMonomial{S<:Statistics}
    factors::Vector{OccupationAtom{S}}

    function OccupationMonomial{S}(
        factors::Vector{OccupationAtom{S}}, ::Val{:raw}
    ) where {S<:Statistics}
        return new{S}(factors)
    end
end

function OccupationMonomial{S}() where {S<:Statistics}
    return OccupationMonomial{S}(OccupationAtom{S}[], Val(:raw))
end

function OccupationMonomial(
    factors::AbstractVector{OccupationAtom{S}}
) where {S<:Statistics}
    canonical = collect(factors)
    sort!(canonical)
    return OccupationMonomial{S}(canonical, Val(:raw))
end

Base.length(monomial::OccupationMonomial) = length(monomial.factors)
Base.isempty(monomial::OccupationMonomial) = isempty(monomial.factors)
Base.iterate(monomial::OccupationMonomial) = iterate(monomial.factors)
Base.iterate(monomial::OccupationMonomial, state) = iterate(monomial.factors, state)
Base.eltype(::Type{OccupationMonomial{S}}) where {S} = OccupationAtom{S}
Base.IteratorSize(::Type{<:OccupationMonomial}) = Base.HasLength()
function Base.isequal(a::OccupationMonomial{S}, b::OccupationMonomial{S}) where {S}
    return isequal(a.factors, b.factors)
end
Base.:(==)(a::OccupationMonomial{S}, b::OccupationMonomial{S}) where {S} = isequal(a, b)
function Base.hash(monomial::OccupationMonomial, h::UInt)
    return hash(monomial.factors, hash(OccupationMonomial, h))
end
function Base.isless(a::OccupationMonomial{S}, b::OccupationMonomial{S}) where {S}
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        ai = a.factors[i]
        bi = b.factors[i]
        isequal(ai, bi) || return isless(ai, bi)
    end
    return length(a) < length(b)
end
function Base.:*(a::OccupationMonomial{S}, b::OccupationMonomial{S}) where {S}
    return OccupationMonomial(vcat(a.factors, b.factors))
end

"""Canonical finite polynomial in exact occupation factors."""
struct OccupationPolynomial{C<:Number,S<:Statistics}
    terms::Vector{Pair{OccupationMonomial{S},C}}

    function OccupationPolynomial{C,S}(
        terms::Vector{Pair{OccupationMonomial{S},C}}, ::Val{:raw}
    ) where {C<:Number,S<:Statistics}
        return new{C,S}(terms)
    end
end

function OccupationPolynomial{C,S}() where {C<:Number,S<:Statistics}
    return OccupationPolynomial{C,S}(Pair{OccupationMonomial{S},C}[], Val(:raw))
end

function _canonical_occupation_polynomial(
    terms::Vector{Pair{OccupationMonomial{S},C}}
) where {C<:Number,S<:Statistics}
    isempty(terms) && return OccupationPolynomial{C,S}()
    sort!(terms; by=first)
    out = Pair{OccupationMonomial{S},C}[]
    sizehint!(out, length(terms))
    for (monomial, coefficient) in terms
        iszero(coefficient) && continue
        if !isempty(out) && isequal(first(out[end]), monomial)
            combined = _simplify(last(out[end]) + coefficient)
            pop!(out)
            iszero(combined) || push!(out, monomial => combined)
        else
            push!(out, monomial => coefficient)
        end
    end
    return OccupationPolynomial{C,S}(out, Val(:raw))
end

function OccupationPolynomial{C,S}(
    terms::Vector{Pair{OccupationMonomial{S},C}}
) where {C<:Number,S<:Statistics}
    return _canonical_occupation_polynomial(copy(terms))
end

function OccupationPolynomial(
    atom::OccupationAtom{S}, coefficient::C
) where {C<:Number,S<:Statistics}
    monomial = OccupationMonomial(OccupationAtom{S}[atom])
    return OccupationPolynomial{C,S}([monomial => coefficient])
end

Base.length(polynomial::OccupationPolynomial) = length(polynomial.terms)
Base.isempty(polynomial::OccupationPolynomial) = isempty(polynomial.terms)
Base.iszero(polynomial::OccupationPolynomial) = isempty(polynomial.terms)
Base.iterate(polynomial::OccupationPolynomial) = iterate(polynomial.terms)
Base.iterate(polynomial::OccupationPolynomial, state) = iterate(polynomial.terms, state)
Base.eltype(::Type{OccupationPolynomial{C,S}}) where {C,S} = Pair{OccupationMonomial{S},C}
Base.IteratorSize(::Type{<:OccupationPolynomial}) = Base.HasLength()
function Base.isequal(
    a::OccupationPolynomial{C,S}, b::OccupationPolynomial{C,S}
) where {C,S}
    return isequal(a.terms, b.terms)
end
Base.:(==)(a::OccupationPolynomial, b::OccupationPolynomial) = isequal(a, b)
function Base.hash(polynomial::OccupationPolynomial, h::UInt)
    return hash(polynomial.terms, hash(OccupationPolynomial, h))
end

function Base.zero(::Type{OccupationPolynomial{C,S}}) where {C<:Number,S<:Statistics}
    return OccupationPolynomial{C,S}()
end
Base.zero(polynomial::OccupationPolynomial{C,S}) where {C,S} = zero(typeof(polynomial))
function Base.one(::Type{OccupationPolynomial{C,S}}) where {C<:Number,S<:Statistics}
    monomial = OccupationMonomial{S}()
    return OccupationPolynomial{C,S}([monomial => one(C)])
end
Base.one(polynomial::OccupationPolynomial{C,S}) where {C,S} = one(typeof(polynomial))

function Base.:+(
    a::OccupationPolynomial{C1,S}, b::OccupationPolynomial{C2,S}
) where {C1<:Number,C2<:Number,S<:Statistics}
    D = promote_type(C1, C2)
    terms = Pair{OccupationMonomial{S},D}[]
    sizehint!(terms, length(a) + length(b))
    for (monomial, coefficient) in a
        push!(terms, monomial => convert(D, coefficient))
    end
    for (monomial, coefficient) in b
        push!(terms, monomial => convert(D, coefficient))
    end
    return _canonical_occupation_polynomial(terms)
end

function Base.:-(polynomial::OccupationPolynomial{C,S}) where {C<:Number,S<:Statistics}
    return (-one(C)) * polynomial
end
Base.:-(a::OccupationPolynomial, b::OccupationPolynomial) = a + (-b)

function Base.:*(
    a::OccupationPolynomial{C1,S}, b::OccupationPolynomial{C2,S}
) where {C1<:Number,C2<:Number,S<:Statistics}
    D = promote_type(C1, C2)
    terms = Pair{OccupationMonomial{S},D}[]
    sizehint!(terms, length(a) * length(b))
    for (ma, ca) in a, (mb, cb) in b
        push!(terms, ma * mb => convert(D, ca) * convert(D, cb))
    end
    return _canonical_occupation_polynomial(terms)
end

function Base.:*(
    coefficient::D, polynomial::OccupationPolynomial{C,S}
) where {D<:Number,C<:Number,S<:Statistics}
    P = promote_type(C, D)
    terms = Pair{OccupationMonomial{S},P}[
        monomial => convert(P, coefficient) * convert(P, value) for
        (monomial, value) in polynomial
    ]
    return _canonical_occupation_polynomial(terms)
end
Base.:*(polynomial::OccupationPolynomial, coefficient::Number) = coefficient * polynomial

"""Statistics sign `σ` in the unified relation `F = 1 + 2σ n`."""
occupation_statistics_sign(::Type{Boson}) = Int8(1)
occupation_statistics_sign(::Type{Fermion}) = Int8(-1)

"""Exact collision normalization `σ/2` relating `I_package` to `C_n`."""
function occupation_collision_factor(::Type{S}) where {S<:Statistics}
    return Int(occupation_statistics_sign(S)) // 2
end

"""Substitute every statistical atom according to `F = 1 + 2σ n`."""
function occupation_substitute(
    polynomial::StatisticalPolynomial{C,S}
) where {C<:Number,S<:Statistics}
    D = promote_type(C, Rational{Int})
    out = zero(OccupationPolynomial{D,S})
    _, slope = statistical_occupation_coefficients(S)
    for (monomial, coefficient) in polynomial
        partial = one(OccupationPolynomial{D,S})
        for atom in monomial
            occupation_atom = OccupationAtom{S}(atom.family, atom.momentum)
            affine =
                one(OccupationPolynomial{D,S}) +
                OccupationPolynomial(occupation_atom, convert(D, slope))
            partial = partial * affine
        end
        out = out + convert(D, coefficient) * partial
    end
    return out
end

"""
Convert a package collision statistical polynomial to the occupation convention.

This performs both the atomwise substitution `F = 1 + 2σ n` and the complete collision
normalization `C_n = σ I_package / 2`.
"""
function occupation_collision_polynomial(
    polynomial::StatisticalPolynomial{C,S}
) where {C<:Number,S<:Statistics}
    expanded = occupation_substitute(polynomial)
    D = promote_type(C, Rational{Int})
    factor = convert(D, occupation_collision_factor(S))
    return factor * expanded
end
