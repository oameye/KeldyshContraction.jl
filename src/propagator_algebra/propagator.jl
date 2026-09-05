#########################
#        Momentum
#########################

struct Momentum
    index::Int8
end
function Momentum(index::Number)
    index < 0 && error("Momentum index must be non-negative, got $index")
    return Momentum(convert(Int8, index))
end
Base.isequal(m1::Momentum, m2::Momentum) = m1.index == m2.index
Base.hash(m::Momentum, h::UInt) = hash(Momenta, hash(m.index, h))

struct Momenta
    prefactors::Vector{Int}
    momenta::Vector{Momentum}
    function Momenta(prefactors::Vector{Int}, momenta::Vector{Momentum}, ::Val{:raw})
        return new(prefactors, momenta)
    end
end
function Momenta(prefactors::Vector{Int}, momenta::Vector{Momentum})
    @assert length(prefactors) == length(momenta) "Length of prefactors and momenta must match"
    idxs = findall(x -> !iszero(x), prefactors)
    return Momenta(prefactors[idxs], momenta[idxs], Val(:raw))
end
Momenta(idx::Int) = Momenta([1], [Momentum(idx)])
Momenta() = Momenta(Int[], Momentum[])

function Base.isequal(m1::Momenta, m2::Momenta)
    return isequal(m1.prefactors, m2.prefactors) && isequal(m1.momenta, m2.momenta)
end
Base.hash(m::Momenta, h::UInt) = hash(Momenta, hash(m.momenta, hash(m.prefactors, h)))

#########################
#      Contraction
#########################

"""Ordered pair of fields forming a two-point contraction."""
struct Contraction{S<:Statistics}
    out::Field{S}
    in::Field{S}
end

Contraction(fields::Tuple{Field{S},Field{S}}) where {S<:Statistics} = Contraction(fields...)
fields(c::Contraction) = (c.out, c.in)

Base.length(::Contraction) = 2
Base.firstindex(::Contraction) = 1
Base.lastindex(::Contraction) = 2
function Base.getindex(c::Contraction, i::Int)
    i == 1 && return c.out
    i == 2 && return c.in
    return throw(BoundsError(c, i))
end
Base.iterate(c::Contraction) = (c.out, 2)
Base.iterate(c::Contraction, state::Int) = state == 2 ? (c.in, 3) : nothing
Base.eltype(::Type{Contraction{S}}) where {S<:Statistics} = Field{S}
Base.IteratorSize(::Type{<:Contraction}) = Base.HasLength()
Base.broadcastable(c::Contraction) = fields(c)
Base.map(f, c::Contraction) = Contraction(f(c.out), f(c.in))
Base.Tuple(c::Contraction) = fields(c)

function Base.isequal(a::Contraction{S}, b::Contraction{S}) where {S<:Statistics}
    return isequal(a.out, b.out) && isequal(a.in, b.in)
end
Base.:(==)(a::Contraction{S}, b::Contraction{S}) where {S<:Statistics} = isequal(a, b)
Base.hash(c::Contraction, h::UInt) = hash(Contraction, hash(c.in, hash(c.out, h)))

function Base.convert(
    ::Type{Contraction{S}}, fields::Tuple{Field{S},Field{S}}
) where {S<:Statistics}
    return Contraction(fields)
end
function Base.convert(
    ::Type{Contraction}, fields::Tuple{Field{S},Field{S}}
) where {S<:Statistics}
    return Contraction(fields)
end

function same_field_family(a::Field{S}, b::Field{S}) where {S<:Statistics}
    return isequal(field_family(a), field_family(b))
end
function contraction_compatible(::Type{Boson}, out::Field{Boson}, in::Field{Boson})
    return same_field_family(out, in)
end

#########################
#         Edge
#########################

"""Bosonic propagator type in the retarded-advanced-Keldysh basis."""
@enumx PropagatorType begin
    Keldysh
    Advanced
    Retarded
    Spectral
end

struct Edge
    out::Field{Boson}
    in::Field{Boson}
    edgetype::PropagatorType.T
    momenta::Momenta
end

function Edge(out::Field{Boson}, in::Field{Boson}, edgetype::PropagatorType.T)
    return Edge(out, in, edgetype, Momenta())
end
Edge(edge::Edge, momenta::Momenta) = Edge(edge.out, edge.in, edge.edgetype, momenta)

function Edge(contraction::Contraction{Boson})
    _out, _in = contraction
    propagator_checks(_out, _in)
    return Edge(_out, _in, propagator_type(_out, _in))
end
Edge(fields::Tuple{Field{Boson},Field{Boson}}) = Edge(Contraction(fields))
Edge(out::Field{Boson}, in::Field{Boson}) = Edge(Contraction(out, in))

momenta(e::Edge) = e.momenta
has_momenta(edge::Edge) = !isempty(edge.momenta.prefactors)

function Base.isequal(e1::Edge, e2::Edge)
    return isequal(e1.out, e2.out) &&
           isequal(e1.in, e2.in) &&
           isequal(e1.edgetype, e2.edgetype)
end
Base.hash(q::Edge, h::UInt) = hash(Edge, hash(q.in, hash(q.edgetype, hash(q.out, h))))

"""Check the rules for a physical bosonic propagator."""
function propagator_checks(out::Field{Boson}, in::Field{Boson})::Nothing
    @assert is_barred(in) "The incoming field must be barred"
    @assert is_unbarred(out) "The outgoing field must be unbarred"
    @assert contraction_compatible(Boson, out, in) "Contracted fields must belong to the same field family"

    v = Contraction(out, in)
    ps = position.(v)
    @assert !is_in(first(ps)) "The outgoing field cannot be at In()"
    @assert !is_out(last(ps)) "The incoming field cannot be at Out()"
    @assert !(has_in(ps) && has_out(ps)) "Cannot contract In() directly with Out()"
    @assert !is_qq_contraction(v) "The quantum-quantum propagator is zero"
    return nothing
end

"""Determine the propagator type from the Keldysh indices of two fields."""
function propagator_type(out::Field{Boson}, in::Field{Boson})::PropagatorType.T
    contours = Int.(keldysh_index.((out, in)))
    diff_contour = first(-(contours...))
    if iszero(diff_contour)
        return PropagatorType.Keldysh
    elseif isone(diff_contour)
        return PropagatorType.Retarded
    else
        return PropagatorType.Advanced
    end
end
function propagator_type(::Type{Boson}, out::Field{Boson}, in::Field{Boson})
    return propagator_type(out, in)
end

propagator_type(e::Edge) = e.edgetype
is_advanced(x::PropagatorType.T) = Int(x) == Int(PropagatorType.Advanced)
is_retarded(x::PropagatorType.T) = Int(x) == Int(PropagatorType.Retarded)
is_keldysh(x::PropagatorType.T) = Int(x) == Int(PropagatorType.Keldysh)
is_spectral(x::PropagatorType.T) = Int(x) == Int(PropagatorType.Spectral)
is_advanced(x::Edge) = is_advanced(propagator_type(x))
is_retarded(x::Edge) = is_retarded(propagator_type(x))
is_keldysh(x::Edge) = is_keldysh(propagator_type(x))
is_spectral(x::Edge) = is_spectral(propagator_type(x))
is_advanced(x::Contraction) = is_advanced(propagator_type(x...))
is_retarded(x::Contraction) = is_retarded(propagator_type(x...))
is_keldysh(x::Contraction) = is_keldysh(propagator_type(x...))

make_spectral(edge::Edge) = Edge(edge.out, edge.in, PropagatorType.Spectral, edge.momenta)
function make_retarded(edge::Edge)
    return Edge(
        bar(edge.in)(position(edge.out)),
        bar(edge.out)(position(edge.in)),
        PropagatorType.Retarded,
        edge.momenta,
    )
end
function make_advanced(edge::Edge)
    return Edge(
        bar(edge.in)(position(edge.out)),
        bar(edge.out)(position(edge.in)),
        PropagatorType.Advanced,
        edge.momenta,
    )
end

fields(e::Edge) = (e.out, e.in)
regularisations(p::Edge) = regularisation.(fields(p))
regularisations(p::Contraction) = regularisation.(p)

function set_reg_to_zero(p::Edge)
    new_fields = map(set_reg_to_zero, fields(p))
    return Edge(new_fields..., p.edgetype, p.momenta)
end
contours(p::Edge) = keldysh_index.(fields(p))

"""Adjoint of a two-point contraction or edge."""
function Base.adjoint(c::Contraction{S}) where {S<:Statistics}
    return Contraction(bar(c.in(position(c.out))), bar(c.out(position(c.in))))
end
Base.adjoint(e::Edge) = Edge(adjoint(Contraction(e.out, e.in)))

"""Reverse a contraction and exchange external coordinates."""
function reverse_contraction(c::Contraction)
    _out, _in = c
    out′, in′ = adjoint(c)
    return Contraction(out′(swap_in_out(position(_in))), in′(swap_in_out(position(_out))))
end
function reverse_edge(e::Edge)
    return Edge(Edge(reverse_contraction(Contraction(e.out, e.in))), momenta(e))
end

#########################
#       Position
#########################

is_bulk(p::Edge) = all(is_bulk, fields(p))
is_bulk(qs::Contraction) = all(is_bulk, qs)
is_in(qs::Contraction) = any(is_in, qs)
is_out(qs::Contraction) = any(is_out, qs)

positions(p::Edge) = position.(fields(p))
positions(p::Contraction) = position.(p)
integer_positions(p::Contraction) = index.(positions(p))
integer_positions(p::Edge) = index.(positions(p))

function same_position(p::Contraction)
    ps = positions(p)
    return isequal(ps...)
end

function position_category(p::Edge)::Symbol
    ps = positions(p)
    if count(is_in, ps) == 1
        return :in
    elseif count(is_out, ps) == 1
        return :out
    elseif all(is_bulk, ps)
        return :bulk
    else
        throw(ArgumentError("Not a valid propagator."))
    end
end

function direction(edge::Edge)
    ps = integer_positions(edge)
    return ps[1] < ps[2]
end
