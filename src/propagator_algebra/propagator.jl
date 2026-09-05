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
#         Edge
#########################

# Transitional tuple representation retained until #241 replaces it with Contraction{S}.
const Contraction = Tuple{Field{Boson},Field{Boson}}

"""
Type of propagator formed by contracting two fields, labelled by their Keldysh indices:
the unbarred field supplies `x` and the barred field `y` in an x-y contour pair.
- `Keldysh`: a Classical-Classical contour;
- `Advanced`: a Quantum-Classical contour;
- `Retarded`: a Classical-Quantum contour;
- `Spectral`: the retarded-minus-advanced combination.

The Quantum-Quantum propagator is always zero, so it has no label here.
"""
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

function Edge(tt::Contraction)
    _out, _in = tt
    propagator_checks(_out, _in)
    return Edge(_out, _in, propagator_type(_out, _in))
end
Edge(out::Field{Boson}, in::Field{Boson}) = Edge((out, in))

momenta(e::Edge) = e.momenta
has_momenta(edge::Edge) = !isempty(edge.momenta.prefactors)

function Base.isequal(e1::Edge, e2::Edge)
    return isequal(e1.out, e2.out) &&
           isequal(e1.in, e2.in) &&
           isequal(e1.edgetype, e2.edgetype)
end
Base.hash(q::Edge, h::UInt) = hash(Edge, hash(q.in, hash(q.edgetype, hash(q.out, h))))

"Collect and check the rules for a physical bosonic propagator."
function propagator_checks(out::Field{Boson}, in::Field{Boson})::Nothing
    @assert is_barred(in) "The incoming field must be barred"
    @assert is_unbarred(out) "The outgoing field must be unbarred"

    v = (out, in)
    ps = position.(v)
    @assert !is_in(first(ps)) "The outgoing field cannot be at In()"
    @assert !is_out(last(ps)) "The incoming field cannot be at Out()"
    @assert !(has_in(ps) && has_out(ps)) "Cannot contract In() directly with Out()"
    @assert !is_qq_contraction(v) "The quantum-quantum propagator is zero"
    return nothing
end

"""
$(DocStringExtensions.SIGNATURES)

Determine the [`PropagatorType`](@ref) in the Retarded-Advanced-Keldysh basis from the
Keldysh indices of the outgoing and incoming fields.
"""
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

"""
Adjoint of a propagator/contraction is meaningful at the two-point-object level: exchange
endpoints and toggle their path-integral orientation. Coordinates are kept in place here;
`reverse_contraction` additionally reverses the coordinates for diagram adjoints.
"""
function Base.adjoint(c::Contraction)
    return (bar(c[2](position(c[1]))), bar(c[1](position(c[2]))))
end
Base.adjoint(c::Edge) = Edge(adjoint(fields(c)))

"""
    reverse_contraction(c::Contraction)

Reverse a contraction while exchanging `In()` and `Out()` coordinates. This is the
coordinate-reversing operation used when taking the adjoint of a whole diagram.
"""
function reverse_contraction(c::Contraction)
    _out, _in = c
    out′, in′ = adjoint(c)
    return (out′(swap_in_out(position(_in))), in′(swap_in_out(position(_out))))
end
reverse_edge(e::Edge) = Edge(Edge(reverse_contraction(fields(e))), momenta(e))

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
