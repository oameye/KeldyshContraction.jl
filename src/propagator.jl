const Contraction = Tuple{<:Destroy,<:Create}

"""
    PropagatorType `Keldysh`, `Advanced`, `Retarded`

The type of propagator taken of two fields with the x-y Contour where `x` is the Contour of the [`Destroy`](@ref) field and `y` the contour of the [`Create`](@ref) field.
- `Keldysh` propagator of a Classical-Classical contour
- `Advanced` propagator of a Quantum-Classical contour
- `Retarded` propagator of a Classical-Quantum contour

The Quantum-Quantum propagator should always be zero.
"""
@enum PropagatorType begin
    Keldysh
    Advanced
    Retarded
end

struct Edge{P2,P1}
    out::Destroy{KeldyshContour,P2,Regularisation}
    in::Create{KeldyshContour,P1,Regularisation}
    edgetype::PropagatorType
    function Edge(tt::Contraction)
        _out, _in = tt[1], tt[2]
        propagator_checks(_out, _in)

        type = propagator_type(_out, _in)
        P2 = position(_out)
        P1 = position(_in)
        return new{typeof(P2),typeof(P1)}(_out, _in, type)
    end
    Edge(out::QSym, in::QSym) = Edge((out, in))
end

"Collect and checks the rules for a physical propagator"
function propagator_checks(out::QSym, in::QSym)::Nothing
    @assert isa(in, Create) "The `in` field must be a Create operator"
    @assert isa(out, Destroy) "The `out` field must be a Destroy operator"
    v = (out, in)

    positions = position.(v)
    @assert !(first(positions) isa In) "The outgoing field can't be the In<:Position` coordinate"
    @assert !(last(positions) isa Out) "The incoming field can't be the Out<:Position` coordinate"
    in_out = (In() ∈ positions ? !(Out() ∈ positions) : true)
    @assert in_out "Can't make a propagator with `In<:Position` and `Out<:Position` coordinate"
    contours = Int.(contour.(v))
    @assert !is_qq_contraction(v) "The quantum-quantum progator is zero"
    return nothing
end

"""
$(TYPEDSIGNATURES)

Determine the type of the propagator in the Retarded-Advance-Keldysh ([`PropagatorType`](@ref)) based on the contour of the output and input quantum field.
"""
function propagator_type(out::QSym, in::QSym)::PropagatorType
    contours = Int.(contour.((out, in)))
    diff_contour = first(-(contours...))
    if iszero(diff_contour)
        return Keldysh
    elseif isone(diff_contour)
        return Retarded
    else
        return Advanced
    end
end
propagator_type(e::Edge) = e.edgetype
is_advanced(x::PropagatorType) = Int(x) == Int(Advanced)
is_retarded(x::PropagatorType) = Int(x) == Int(Retarded)
is_keldysh(x::PropagatorType) = Int(x) == Int(Keldysh)
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

isbulk(p::Edge) = all(isbulk.(fields(p)))
isbulk(qs::Contraction) = all(isbulk.(qs))
is_in(qs::Contraction) = any(is_in.(qs))
is_out(qs::Contraction) = any(is_out.(qs))

function positions(p::Edge)
    return position.(fields(p))
end
function positions(p::Contraction)
    return position.(p)
end

function integer_positions(p::Contraction)::NTuple{2,Int}
    return Int.(positions(p))
end
function integer_positions(p::Edge)::NTuple{2,Int}
    return Int.(positions(p))
end
function same_position(p::Contraction)
    _positions = positions(p)
    return isequal(_positions...)
end

function position(p::Edge)
    _positions = positions(p)
    if length(findall(x -> x isa In, _positions)) == 1
        return In()
    elseif length(findall(x -> x isa Out, _positions)) == 1
        return Out()
    elseif all(isbulk, _positions)
        idxs = map(x -> getproperty(x, :index), _positions)
        if isequal(idxs...)
            return _positions[1]
        else
            return minimum(_positions)
            # ^  TODO what to do if both are bulk with different index?
        end
    else
        throw(ArgumentError("Not a valid propagator."))
    end
end # TODO move to Hopcroft-Ullman pairing
