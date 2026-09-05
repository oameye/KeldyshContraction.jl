#########################
#       Field{S}
#########################

"""
    Field{S} <: QSym

Fundamental path-integral field with statistics `S`. Statistics is encoded in the type;
orientation, Keldysh index, position, regularisation, and internal indices are concrete
runtime values.
"""
struct Field{S<:Statistics} <: QSym
    name::Symbol
    orientation::Orientation.T
    keldysh::KeldyshIndex.T
    position::Position
    regularisation::Regularisation.T
    indices::FieldIndices
end

function Field{S}(
    name::Symbol,
    keldysh::KeldyshIndex.T,
    orientation::Orientation.T=Orientation.Unbarred,
    reg::Regularisation.T=Regularisation.Zero,
    pos::Position=Bulk(),
    indices::FieldIndices=FieldIndices(),
) where {S<:Statistics}
    return Field{S}(name, orientation, keldysh, pos, reg, indices)
end

statistics(::Field{S}) where {S<:Statistics} = S
orientation(f::Field) = f.orientation
keldysh_index(f::Field) = f.keldysh
contour(f::Field) = keldysh_index(f)
regularisation(f::Field) = f.regularisation
position(f::Field) = f.position
field_indices(f::Field) = f.indices
index(f::Field) = index(position(f))

is_unbarred(f::Field) = orientation(f) === Orientation.Unbarred
is_barred(f::Field) = orientation(f) === Orientation.Barred
is_annihilation(f::Field) = is_unbarred(f)
is_creation(f::Field) = is_barred(f)
ladder(f::Field) = Int(is_barred(f))

is_quantum(f::Field{Boson}) = keldysh_index(f) === Quantum
is_classical(f::Field{Boson}) = keldysh_index(f) === Classical

is_bulk(f::Field) = is_bulk(position(f))
is_in(f::Field) = is_in(position(f))
is_out(f::Field) = is_out(position(f))

function (f::Field{S})(pos::Position) where {S}
    return Field{S}(
        name(f), orientation(f), keldysh_index(f), pos, regularisation(f), field_indices(f)
    )
end
function (f::Field{S})(reg::Regularisation.T) where {S}
    return Field{S}(
        name(f), orientation(f), keldysh_index(f), position(f), reg, field_indices(f)
    )
end

function set_reg_to_zero(f::Field{S}) where {S}
    return Field{S}(
        name(f),
        orientation(f),
        keldysh_index(f),
        position(f),
        Regularisation.Zero,
        field_indices(f),
    )
end
function set_reg_to_zero!(v::Vector{Field{S}}) where {S}
    for i in eachindex(v)
        v[i] = set_reg_to_zero(v[i])
    end
    return v
end

contour_integers(v::Vector{Field{S}}) where {S} = Int[Int(keldysh_index(x)) for x in v]

"""
    bar(field)

Toggle the independent path-integral orientation of a field. `bar` is deliberately
separate from `adjoint`: a barred integration variable is not represented as a field-level
Hermitian adjoint operation.
"""
function bar(f::Field{S}) where {S}
    o = is_unbarred(f) ? Orientation.Barred : Orientation.Unbarred
    return Field{S}(
        name(f), o, keldysh_index(f), position(f), regularisation(f), field_indices(f)
    )
end

function Base.isequal(a::Field{S}, b::Field{S}) where {S}
    return isequal(name(a), name(b)) &&
           isequal(orientation(a), orientation(b)) &&
           isequal(keldysh_index(a), keldysh_index(b)) &&
           isequal(position(a), position(b)) &&
           isequal(regularisation(a), regularisation(b)) &&
           isequal(field_indices(a), field_indices(b))
end
Base.:(==)(a::Field{S}, b::Field{S}) where {S} = isequal(a, b)

function Base.isless(a::Field{S}, b::Field{S}) where {S}
    oa, ob = Int(orientation(a)), Int(orientation(b))
    oa == ob || return oa < ob
    pa, pb = index(a), index(b)
    pa == pb || return pa < pb
    ka, kb = Int(keldysh_index(a)), Int(keldysh_index(b))
    ka == kb || return ka < kb
    name(a) == name(b) || return isless(name(a), name(b))
    ra, rb = Int(regularisation(a)), Int(regularisation(b))
    ra == rb || return ra < rb
    return isless(field_indices(a), field_indices(b))
end

"""Statistics-dependent sign for exchanging two fields during canonical ordering."""
exchange_sign(::Type{Boson}, ::Field{Boson}, ::Field{Boson}) = Int8(1)

"""
Canonicalize a concrete field vector in place and return the accumulated exchange sign.
The insertion-sort implementation is intentional: future fermionic ordering can account
for every actual exchange without changing the storage or ordering algorithm.
"""
function canonicalize_fields!(args::Vector{Field{S}}) where {S<:Statistics}
    sign = Int8(1)
    for i in 2:length(args)
        j = i
        while j > 1 && isless(args[j], args[j - 1])
            sign *= exchange_sign(S, args[j - 1], args[j])
            args[j - 1], args[j] = args[j], args[j - 1]
            j -= 1
        end
    end
    return sign
end

Base.one(f::Field{S}) where {S<:Statistics} = QMul{Int,S}(1, Field{S}[])
Base.zero(f::Field{S}) where {S<:Statistics} = QMul{Int,S}(0, Field{S}[])
Base.one(::Type{Field{S}}) where {S<:Statistics} = QMul{Int,S}(1, Field{S}[])
Base.zero(::Type{Field{S}}) where {S<:Statistics} = QMul{Int,S}(0, Field{S}[])

"""
    @qfields

Construct path-integral fields. The statistics type is supplied before the semantic
Keldysh label:

```julia
@qfields c::Boson(Classical) q::Boson(Quantum)
bar(c)
```
"""
macro qfields(qs...)
    defs = map(qs) do q
        nf = _name_field(q)
        fname, field_expr = nf.name, nf.field_expr
        statistics_expr = field_expr.args[1]
        field_args = field_expr.args[2:end]

        construction = :(Field{$(esc(statistics_expr))}(
            $(QuoteNode(fname)), $(map(esc, field_args)...), Orientation.Unbarred
        ))
        return :($(esc(fname)) = $construction)
    end

    names = map(q -> esc(_name_field(q).name), qs)
    return Expr(:block, defs..., :(tuple($(names...))))
end

function _name_field(expr)
    @assert expr isa Expr && expr.head == :(::) "Expected expression of form name::Statistics(args...)"
    name = expr.args[1]
    @assert name isa Symbol "Left side of :: must be a symbol"

    field_expr = expr.args[2]
    @assert field_expr isa Expr && field_expr.head == :call "Right side of :: must be a field constructor call"

    return (name=name, field_expr=field_expr)
end
