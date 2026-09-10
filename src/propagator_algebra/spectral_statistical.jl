"""Causal/spectral character of one line in the post-Wigner kinetic IR."""
@enum KineticLineKind::UInt8 begin
    KineticRetarded = 0
    KineticAdvanced = 1
    KineticSpectral = 2
end

"""Whether a spectral line carries the statistical distribution factor `F`."""
@enum StatisticalWeight::UInt8 begin
    NoStatisticalWeight = 0
    DistributionWeight = 1
end

"""
One physical propagator line after exact R/A/K-to-spectral/statistical lowering.

The line no longer stores a Keldysh matrix index. It retains the physical field family,
endpoint positions, the canonical relative equal-time regularisation shift, exact routed
momentum, causal/spectral kind, and whether the line carries the statistical factor `F`.
"""
struct KineticLine{S<:Statistics}
    family::FieldFamily{S}
    out_position::Position
    in_position::Position
    regularisation_shift::Int8
    momentum::LinearMomentum
    kind::KineticLineKind
    statistical::StatisticalWeight
end

statistics(::KineticLine{S}) where {S} = S
kinetic_line_kind(line::KineticLine) = line.kind
statistical_weight(line::KineticLine) = line.statistical
regularisation_shift(line::KineticLine) = line.regularisation_shift
momentum(line::KineticLine) = line.momentum

function Base.isequal(a::KineticLine{S}, b::KineticLine{S}) where {S<:Statistics}
    return isequal(a.family, b.family) &&
           isequal(a.out_position, b.out_position) &&
           isequal(a.in_position, b.in_position) &&
           a.regularisation_shift == b.regularisation_shift &&
           isequal(a.momentum, b.momentum) &&
           a.kind === b.kind &&
           a.statistical === b.statistical
end
Base.:(==)(a::KineticLine{S}, b::KineticLine{S}) where {S<:Statistics} = isequal(a, b)
function Base.hash(line::KineticLine, h::UInt)
    h = hash(KineticLine, h)
    h = hash(line.family, h)
    h = hash(line.out_position, h)
    h = hash(line.in_position, h)
    h = hash(line.regularisation_shift, h)
    h = hash(line.momentum, h)
    h = hash(line.kind, h)
    return hash(line.statistical, h)
end

function Base.isless(a::KineticLine{S}, b::KineticLine{S}) where {S<:Statistics}
    isequal(a.family, b.family) || return isless(a.family, b.family)
    isequal(a.out_position, b.out_position) || return isless(a.out_position, b.out_position)
    isequal(a.in_position, b.in_position) || return isless(a.in_position, b.in_position)
    a.regularisation_shift == b.regularisation_shift ||
        return a.regularisation_shift < b.regularisation_shift
    isequal(a.momentum, b.momentum) ||
        return _linear_momentum_isless(a.momentum, b.momentum)
    a.kind === b.kind || return Int(a.kind) < Int(b.kind)
    return Int(a.statistical) < Int(b.statistical)
end

@inline function _kinetic_regularisation_shift(edge::Edge)::Int8
    out_field, in_field = fields(edge)
    same_position(Contraction(out_field, in_field)) || return Int8(0)
    propagator = propagator_type(edge)
    (is_retarded(propagator) || is_advanced(propagator)) || return Int8(0)
    return convert(Int8, subtraction((regularisation(out_field), regularisation(in_field))))
end

@inline function _kinetic_line(
    edge::Edge{S},
    routed_momentum::LinearMomentum,
    kind::KineticLineKind,
    statistical::StatisticalWeight,
) where {S<:Statistics}
    out_field, in_field = fields(edge)
    return KineticLine{S}(
        field_family(out_field),
        position(out_field),
        position(in_field),
        _kinetic_regularisation_shift(edge),
        routed_momentum,
        kind,
        statistical,
    )
end

"""Canonical commutative product of physical kinetic lines."""
struct KineticMonomial{S<:Statistics,E}
    lines::FixedVector{E,KineticLine{S}}
end

function KineticMonomial(lines::Vector{KineticLine{S}}, ::Val{E}) where {S<:Statistics,E}
    length(lines) == E ||
        throw(ArgumentError("kinetic monomial line count does not match Val{$E}"))
    canonical = copy(lines)
    sort!(canonical)
    return KineticMonomial{S,E}(FixedVector{E,KineticLine{S}}(canonical))
end

kinetic_lines(monomial::KineticMonomial) = monomial.lines
Base.length(monomial::KineticMonomial) = length(monomial.lines)
Base.iterate(monomial::KineticMonomial) = iterate(monomial.lines)
Base.iterate(monomial::KineticMonomial, state) = iterate(monomial.lines, state)
Base.eltype(::Type{KineticMonomial{S,E}}) where {S,E} = KineticLine{S}
Base.IteratorSize(::Type{<:KineticMonomial}) = Base.HasLength()
function Base.isequal(a::KineticMonomial{S,E}, b::KineticMonomial{S,E}) where {S,E}
    return isequal(a.lines, b.lines)
end
Base.:(==)(a::KineticMonomial{S,E}, b::KineticMonomial{S,E}) where {S,E} = isequal(a, b)
function Base.hash(monomial::KineticMonomial, h::UInt)
    return hash(KineticMonomial, hash(monomial.lines, h))
end

function Base.:*(a::KineticLine{S}, b::KineticLine{S}) where {S<:Statistics}
    return KineticMonomial(KineticLine{S}[a, b], Val(2))
end

function Base.:*(
    a::KineticMonomial{S,E1}, b::KineticMonomial{S,E2}
) where {S<:Statistics,E1,E2}
    lines = KineticLine{S}[]
    sizehint!(lines, E1 + E2)
    append!(lines, a.lines)
    append!(lines, b.lines)
    return KineticMonomial(lines, Val(E1 + E2))
end

function Base.:*(
    monomial::KineticMonomial{S,E}, line::KineticLine{S}
) where {S<:Statistics,E}
    lines = KineticLine{S}[]
    sizehint!(lines, E + 1)
    append!(lines, monomial.lines)
    push!(lines, line)
    return KineticMonomial(lines, Val(E + 1))
end
function Base.:*(
    line::KineticLine{S}, monomial::KineticMonomial{S,E}
) where {S<:Statistics,E}
    return monomial * line
end

"""
Canonical kinetic term for one routed graph contribution.

The propagator product is stored as a canonical `KineticMonomial`. `basis` and
`external_momentum` preserve the exact meaning of every `LinearMomentum`. The derivative
vertex factor remains a separate exact `MomentumPolynomial`.
"""
struct KineticTerm{S<:Statistics,E1,E2}
    monomial::KineticMonomial{S,E1}
    topology::FixedVector{E2,Int}
    basis::MomentumBasis
    external_momentum::MomentumVariable
    kinematic::MomentumPolynomial{ComplexRationals}
end

statistics(::KineticTerm{S}) where {S} = S
kinetic_lines(term::KineticTerm) = kinetic_lines(term.monomial)
kinematic_factor(term::KineticTerm) = term.kinematic
momentum_basis(term::KineticTerm) = term.basis
external_wigner_momentum(term::KineticTerm) = term.external_momentum

@inline function _kinetic_topology_isequal(
    a::FixedVector{E,Int}, b::FixedVector{E,Int}
) where {E}
    @inbounds for i in 1:E
        a[i] == b[i] || return false
    end
    return true
end

@inline function _kinetic_topology_hash(topology::FixedVector{E,Int}, h::UInt) where {E}
    h = hash(E, h)
    @inbounds for i in 1:E
        h = hash(topology[i], h)
    end
    return h
end

function Base.isequal(
    a::KineticTerm{S,E1,E2}, b::KineticTerm{S,E1,E2}
) where {S<:Statistics,E1,E2}
    return isequal(a.monomial, b.monomial) &&
           _kinetic_topology_isequal(a.topology, b.topology) &&
           isequal(a.basis, b.basis) &&
           isequal(a.external_momentum, b.external_momentum) &&
           isequal(a.kinematic, b.kinematic)
end
Base.:(==)(a::KineticTerm{S,E1,E2}, b::KineticTerm{S,E1,E2}) where {S,E1,E2} = isequal(a, b)
function Base.hash(term::KineticTerm{S,E1,E2}, h::UInt) where {S,E1,E2}
    h = hash(KineticTerm, h)
    h = hash(term.monomial, h)
    h = _kinetic_topology_hash(term.topology, h)
    h = hash(term.basis, h)
    h = hash(term.external_momentum, h)
    return hash(term.kinematic, h)
end

"""Concrete spectral/statistical polynomial at fixed statistics, graph shape, and Wigner order."""
struct KineticExpression{C<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext}
    terms::Dict{KineticTerm{S,E1,E2},C}
    context::Ctx
end

function KineticExpression{C,S,E1,E2,G,Ctx}(
    context::Ctx
) where {C<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext}
    K = KineticTerm{S,E1,E2}
    return KineticExpression{C,S,E1,E2,G,Ctx}(Dict{K,C}(), context)
end

statistics(::KineticExpression{C,S}) where {C,S} = S
gradient_order(::KineticExpression{C,S,E1,E2,G}) where {C,S,E1,E2,G} = Val(G)
wigner_context(expression::KineticExpression) = expression.context
Base.length(expression::KineticExpression) = length(expression.terms)
Base.isempty(expression::KineticExpression) = isempty(expression.terms)
Base.iszero(expression::KineticExpression) = isempty(expression.terms)
Base.iterate(expression::KineticExpression) = iterate(expression.terms)
Base.iterate(expression::KineticExpression, state) = iterate(expression.terms, state)
function Base.eltype(::Type{KineticExpression{C,S,E1,E2,G,Ctx}}) where {C,S,E1,E2,G,Ctx}
    return Pair{KineticTerm{S,E1,E2},C}
end

function Base.isequal(
    a::KineticExpression{C,S,E1,E2,G,Ctx}, b::KineticExpression{C,S,E1,E2,G,Ctx}
) where {C,S,E1,E2,G,Ctx}
    return isequal(a.terms, b.terms) && isequal(a.context, b.context)
end
Base.:(==)(a::KineticExpression, b::KineticExpression) = isequal(a, b)
function Base.hash(expression::KineticExpression, h::UInt)
    return hash(expression.context, hash(expression.terms, hash(KineticExpression, h)))
end

function Base.push!(
    expression::KineticExpression{C,S,E1,E2,G,Ctx},
    term::KineticTerm{S,E1,E2},
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

function Base.:+(
    a::KineticExpression{C1,S,E1,E2,G,Ctx}, b::KineticExpression{C2,S,E1,E2,G,Ctx}
) where {C1<:Number,C2<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext}
    isequal(a.context, b.context) || throw(
        ArgumentError("cannot add kinetic expressions with different Wigner contexts")
    )
    D = promote_type(C1, C2)
    out = KineticExpression{D,S,E1,E2,G,Ctx}(a.context)
    for (term, coefficient) in a
        push!(out, term, convert(D, coefficient))
    end
    for (term, coefficient) in b
        push!(out, term, convert(D, coefficient))
    end
    return out
end

function Base.:*(
    prefactor::P, expression::KineticExpression{C,S,E1,E2,G,Ctx}
) where {P<:Number,C<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext}
    D = promote_type(P, C)
    out = KineticExpression{D,S,E1,E2,G,Ctx}(expression.context)
    p = convert(D, prefactor)
    for (term, coefficient) in expression
        push!(out, term, p * convert(D, coefficient))
    end
    return out
end
Base.:*(expression::KineticExpression, prefactor::Number) = prefactor * expression
function Base.:-(expression::KineticExpression{C}) where {C<:Number}
    return -one(C) * expression
end
Base.:-(a::KineticExpression, b::KineticExpression) = a + (-b)

kinetic_coefficient_type(::Type{C}) where {C<:Number} = promote_type(C, ComplexRationals)

@inline function _lower_kinetic_line(
    edge::Edge{S}, routed_momentum::LinearMomentum
) where {S<:Statistics}
    propagator = propagator_type(edge)
    if is_keldysh(propagator)
        return _kinetic_line(edge, routed_momentum, KineticSpectral, DistributionWeight),
        Int8(1)
    elseif is_retarded(propagator)
        return _kinetic_line(edge, routed_momentum, KineticRetarded, NoStatisticalWeight),
        Int8(0)
    elseif is_advanced(propagator)
        return _kinetic_line(edge, routed_momentum, KineticAdvanced, NoStatisticalWeight),
        Int8(0)
    elseif is_spectral(propagator)
        return _kinetic_line(edge, routed_momentum, KineticSpectral, NoStatisticalWeight),
        Int8(0)
    end
    return error("unsupported propagator type in kinetic lowering")
end

function _lower_kinetic_term(
    graph::WignerDiagram{S,E1,E2,G,K,Ctx}, contribution::WignerContribution{C}, ::Type{D}
) where {C<:Number,D<:Number,S<:Statistics,E1,E2,G,K,Ctx<:AbstractWignerContext}
    coordinate = coordinate_diagram(graph)
    source_edges = contractions(coordinate)
    routed_momenta = edge_momenta(graph)
    lines = Vector{KineticLine{S}}(undef, E1)
    keldysh_count = Int8(0)
    @inbounds for i in 1:E1
        line, is_statistical = _lower_kinetic_line(source_edges[i], routed_momenta[i])
        lines[i] = line
        keldysh_count += is_statistical
    end

    term = KineticTerm{S,E1,E2}(
        KineticMonomial(lines, Val(E1)),
        topology(coordinate),
        momentum_basis(graph),
        external_wigner_momentum(graph),
        contribution.kinematic,
    )
    phase = convert(D, (-im)^Int(keldysh_count))
    coefficient = _simplify(convert(D, contribution.coefficient) * phase)
    return term, coefficient
end

"""
    kinetic_expression(diagrams::WignerDiagrams)

Lower one complete homogeneous Wigner R/A/K component to the exact spectral/statistical kinetic
IR. The only per-line basis identity applied here is `Gᴷ = -im * F * A`. Retarded and advanced
lines remain retarded and advanced. In particular, this function never approximates the
imaginary part of a multi-line product by multiplying imaginary parts of individual lines.

Only gradient order zero is supported. Higher Wigner gradient orders must retain the noncommuting
Moyal products and therefore require a different kinetic lowering rule.
"""
function kinetic_expression(
    diagrams::WignerDiagrams{C,S,E1,E2,0,Ctx}
) where {C<:Number,S<:Statistics,E1,E2,Ctx<:AbstractWignerContext}
    D = kinetic_coefficient_type(C)
    out = KineticExpression{D,S,E1,E2,0,Ctx}(wigner_context(diagrams))
    for (graph, contributions) in diagrams
        for contribution in contributions
            term, coefficient = _lower_kinetic_term(graph, contribution, D)
            push!(out, term, coefficient)
        end
    end
    return out
end

function kinetic_expression(
    ::WignerDiagrams{C,S,E1,E2,G,Ctx}
) where {C<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext}
    return throw(
        ArgumentError(
            "spectral/statistical kinetic lowering supports only Wigner gradient order zero; got $G",
        ),
    )
end

"""Exact post-Wigner self-energy with semantic R/A/K components in one kinetic IR."""
struct KineticSelfEnergy{C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    keldysh::KineticExpression{C,S,E1,E2,G,Ctx}
    retarded::KineticExpression{C,S,E1,E2,G,Ctx}
    advanced::KineticExpression{C,S,E1,E2,G,Ctx}
    parameter::ParameterMonomial
    target::FieldFamily{S}
    context::Ctx
end

order(::KineticSelfEnergy{C,S,O}) where {C,S,O} = O
statistics(::KineticSelfEnergy{C,S}) where {C,S} = S
parameters(Σ::KineticSelfEnergy) = Σ.parameter
target_family(Σ::KineticSelfEnergy) = Σ.target
gradient_order(::KineticSelfEnergy{C,S,O,E1,E2,G}) where {C,S,O,E1,E2,G} = Val(G)
wigner_context(Σ::KineticSelfEnergy) = Σ.context

function Base.isequal(
    a::KineticSelfEnergy{C,S,O,E1,E2,G,Ctx}, b::KineticSelfEnergy{C,S,O,E1,E2,G,Ctx}
) where {C,S,O,E1,E2,G,Ctx}
    return isequal(a.keldysh, b.keldysh) &&
           isequal(a.retarded, b.retarded) &&
           isequal(a.advanced, b.advanced) &&
           isequal(a.parameter, b.parameter) &&
           isequal(a.target, b.target) &&
           isequal(a.context, b.context)
end
Base.:(==)(a::KineticSelfEnergy, b::KineticSelfEnergy) = isequal(a, b)
function Base.hash(Σ::KineticSelfEnergy, h::UInt)
    return hash((Σ.keldysh, Σ.retarded, Σ.advanced, Σ.parameter, Σ.target, Σ.context), h)
end

function kinetic_expression(
    Σ::WignerSelfEnergy{C,S,O,E1,E2,0,Ctx}
) where {C<:Number,S<:Statistics,O,E1,E2,Ctx<:AbstractWignerContext}
    D = kinetic_coefficient_type(C)
    return KineticSelfEnergy{D,S,O,E1,E2,0,Ctx}(
        kinetic_expression(Σ.keldysh),
        kinetic_expression(Σ.retarded),
        kinetic_expression(Σ.advanced),
        Σ.parameter,
        Σ.target,
        Σ.context,
    )
end

function kinetic_expression(
    ::WignerSelfEnergy{C,S,O,E1,E2,G,Ctx}
) where {C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    return throw(
        ArgumentError(
            "spectral/statistical kinetic lowering supports only Wigner gradient order zero; got $G",
        ),
    )
end

"""Form the complete self-energy discontinuity `Σᴿ - Σᴬ` before any on-shell reduction."""
function retarded_minus_advanced(Σ::KineticSelfEnergy)
    return Σ.retarded - Σ.advanced
end

"""Spectral self-energy in the package convention `A_Σ = im * (Σᴿ - Σᴬ)`."""
function spectral_self_energy(Σ::KineticSelfEnergy)
    return im * retarded_minus_advanced(Σ)
end

"""Affine coefficients `(offset, slope)` in `F = offset + slope*n`."""
statistical_occupation_coefficients(::Type{Boson}) = (Int8(1), Int8(2))
statistical_occupation_coefficients(::Type{Fermion}) = (Int8(1), Int8(-2))

function statistical_from_occupation(::Type{S}, n::N) where {S<:Statistics,N<:Number}
    offset, slope = statistical_occupation_coefficients(S)
    return convert(N, offset) + convert(N, slope) * n
end
