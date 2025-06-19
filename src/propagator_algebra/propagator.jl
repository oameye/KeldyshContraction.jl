#########################
#        Momentum
#########################

struct Momentum
    index::Int8
    function Momentum(index)
        index < 0 && error("Momentum index must be non-negative, got $index")
        return new(index)
    end
end
Base.isequal(m1::Momentum, m2::Momentum) = m1.index == m2.index
struct Momenta
    prefactors::Vector{Int}
    momenta::Vector{Momentum}
    function Momenta(prefactors::Vector{Int}, momenta::Vector{Momentum})
        @assert length(prefactors) == length(momenta) "Length of prefactors and momenta must match"
        idxs = findall(x -> !iszero(x), prefactors)
        return new(prefactors[idxs], momenta[idxs])
    end
end
Momenta(idx::Int) = Momenta([1], [Momentum(idx)])
Momenta() = Momenta(Int[], Momentum[])

function Base.isequal(m1::Momenta, m2::Momenta)
    isequal(m1.prefactors, m2.prefactors) && isequal(m1.momenta, m2.momenta)
end
# SmallCollections.default(::Type{Momenta}) = Momenta(Int[], Momentum[])

#########################
#         Edge
#########################

const Contraction = Tuple{<:Destroy,<:Create}

"""
    PropagatorType `Keldysh`, `Advanced`, `Retarded`

The type of propagator taken of two fields with the x-y Contour where `x` is the Contour of the [`Destroy`](@ref) field and `y` the contour of the [`Create`](@ref) field.
- `Keldysh` propagator of a Classical-Classical contour
- `Advanced` propagator of a Quantum-Classical contour
- `Retarded` propagator of a Classical-Quantum contour

The Quantum-Quantum propagator should always be zero.
"""
@enumx PropagatorType begin
    Keldysh
    Advanced
    Retarded
end

struct Edge
    out::Destroy
    in::Create
    edgetype::PropagatorType.T
    momenta::Momenta

    function Edge(out::Destroy, in::Create, edgetype::PropagatorType.T)
        new(out, in, edgetype, Momenta())
    end
    Edge(edge::Edge, momenta::Momenta) = new(edge.out, edge.in, edge.edgetype, momenta)
end
function Edge(tt::Contraction)
    _out, _in = tt[1], tt[2]
    propagator_checks(_out, _in)

    type = propagator_type(_out, _in)
    P2 = position(_out)
    P1 = position(_in)
    return Edge(_out, _in, type)
end
Edge(out::QSym, in::QSym) = Edge((out, in))

momenta(e::Edge) = e.momenta
has_momenta(edge::Edge) = !isempty(edge.momenta.prefactors)

function Base.isequal(e1::Edge, e2::Edge)
    isequal(e1.out, e2.out) && isequal(e1.in, e2.in) && isequal(e1.edgetype, e2.edgetype)
end
Base.hash(q::Edge, h::UInt) = hash(Edge, hash(q.in, hash(q.edgetype, hash(q.out, h))))

"Collect and checks the rules for a physical propagator"
function propagator_checks(out::QSym, in::QSym)::Nothing
    @assert isa(in, Create) "The `in` field must be a Create operator"
    @assert isa(out, Destroy) "The `out` field must be a Destroy operator"
    v = (out, in)

    positions = position.(v)
    @assert !(is_in(first(positions))) "The outgoing field can't be the In<:Position` coordinate"
    @assert !(is_out(last(positions))) "The incoming field can't be the Out<:Position` coordinate"
    in_out = (has_in(positions) ? !(has_out(positions)) : true)
    @assert in_out "Can't make a propagator with `In<:Position` and `Out<:Position` coordinate"
    contours = Int.(contour.(v))
    @assert !is_qq_contraction(v) "The quantum-quantum progator is zero"
    return nothing
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Determine the type of the propagator in the Retarded-Advance-Keldysh ([`PropagatorType`](@ref)) based on the contour of the output and input quantum field.
"""
function propagator_type(out::QSym, in::QSym)::PropagatorType.T
    contours = Int.(contour.((out, in)))
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
is_advanced(x::Edge) = is_advanced(propagator_type(x))
is_retarded(x::Edge) = is_retarded(propagator_type(x))
is_keldysh(x::Edge) = is_keldysh(propagator_type(x))
is_advanced(x::Contraction) = is_advanced(propagator_type(x...))
is_retarded(x::Contraction) = is_retarded(propagator_type(x...))
is_keldysh(x::Contraction) = is_keldysh(propagator_type(x...))

fields(e::Edge) = (e.out, e.in)
function regularisations(p::Edge)
    return regularisation.(fields(p))
end
function regularisations(p::Contraction)
    return regularisation.(p)
end
contours(p::Edge) = contour.(fields(p))

"""
Gᴷ(x₁, x₂)† = -1*Gᴷ(x₁, x₂)
Gᴬ(x₁, x₂)† = Gᴿ(x₁, x₂)
Gᴿ(x₁, x₂)† = Gᴬ(x₁, x₂)
"""
Base.adjoint(c::Contraction) = adjoint.((c[2](position(c[1])), c[1](position(c[2]))))
Base.adjoint(c::Edge) = Edge(adjoint(fields(c)))

#########################
#       Position
#########################

is_bulk(p::Edge) = all(is_bulk.(fields(p)))
is_bulk(qs::Contraction) = all(is_bulk.(qs))
is_in(qs::Contraction) = any(is_in.(qs))
is_out(qs::Contraction) = any(is_out.(qs))

function positions(p::Edge)
    return position.(fields(p))
end
function positions(p::Contraction)
    return position.(p)
end

function integer_positions(p::Contraction)
    return index.(positions(p))
end
function integer_positions(p::Edge)
    return index.(positions(p))
end
function same_position(p::Contraction)
    _positions = positions(p)
    return isequal(_positions...)
end

function position_category(p::Edge)::Symbol
    _positions = positions(p)
    if length(findall(is_in, _positions)) == 1
        return :in
    elseif length(findall(is_out, _positions)) == 1
        return :out
    elseif all(is_bulk, _positions)
        return :bulk
    else
        throw(ArgumentError("Not a valid propagator."))
    end
end
