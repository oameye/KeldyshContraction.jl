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

"""
$(DocStringExtensions.SIGNATURES)

Construct a field from its Keldysh index, with every other property defaulted. Argument
order here differs from the struct's field order, so prefer keywords past `keldysh`; to
derive one field from another use [`reconstruct`](@ref) rather than restating all six
properties positionally.
"""
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

"""
$(DocStringExtensions.SIGNATURES)

Copy a field, replacing only the named properties and preserving statistics. Every
orientation, position, and regularisation update goes through here, so the six-property
argument clump is written once.
"""
function reconstruct(
    f::Field{S};
    name::Symbol=name(f),
    orientation::Orientation.T=orientation(f),
    keldysh::KeldyshIndex.T=keldysh_index(f),
    position::Position=position(f),
    regularisation::Regularisation.T=regularisation(f),
    indices::FieldIndices=field_indices(f),
) where {S<:Statistics}
    return Field{S}(name, orientation, keldysh, position, regularisation, indices)
end

statistics(::Field{S}) where {S<:Statistics} = S
orientation(f::Field) = f.orientation
keldysh_index(f::Field) = f.keldysh
regularisation(f::Field) = f.regularisation
position(f::Field) = f.position
field_indices(f::Field) = f.indices
index(f::Field) = index(position(f))

is_unbarred(f::Field) = orientation(f) === Orientation.Unbarred
is_barred(f::Field) = orientation(f) === Orientation.Barred
ladder(f::Field) = Int(is_barred(f))

is_quantum(f::Field{Boson}) = keldysh_index(f) === Quantum
is_classical(f::Field{Boson}) = keldysh_index(f) === Classical

is_bulk(f::Field) = is_bulk(position(f))
is_in(f::Field) = is_in(position(f))
is_out(f::Field) = is_out(position(f))

(f::Field)(pos::Position) = reconstruct(f; position=pos)
(f::Field)(reg::Regularisation.T) = reconstruct(f; regularisation=reg)

set_reg_to_zero(f::Field) = reconstruct(f; regularisation=Regularisation.Zero)
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
function bar(f::Field)
    o = is_unbarred(f) ? Orientation.Barred : Orientation.Unbarred
    return reconstruct(f; orientation=o)
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

"""
    exchange_sign(::Type{S}) -> Int8

Sign picked up when two adjacent fields of statistics `S` are exchanged. Bosonic exchange is
sign-free. Adding fermions means adding `exchange_sign(::Type{Fermion}) = Int8(-1)` and
nothing else: canonical ordering, the sign-free fast path, and product construction all read
this one method.
"""
exchange_sign(::Type{Boson}) = Int8(1)

"""
    exchange_sign(::Type{S}, a::Field{S}, b::Field{S}) -> Int8

Sign for exchanging two specific fields. Defers to the statistics-level
[`exchange_sign`](@ref); a statistics needing a pair-dependent rule overrides this method.
"""
function exchange_sign(::Type{S}, ::Field{S}, ::Field{S}) where {S<:Statistics}
    return exchange_sign(S)
end

"""
    is_exchange_sign_free(::Type{S}) -> Bool

Whether every exchange under `S` carries sign `+1`. Derived from [`exchange_sign`](@ref) so
the two can never disagree, and constant-folded at compile time.
"""
is_exchange_sign_free(::Type{S}) where {S<:Statistics} = isone(exchange_sign(S))

"""
    canonicalize_fields!(args::Vector{Field{S}}) -> Int8

Sort a concrete field vector in place and return the accumulated exchange sign.

Sign-free statistics take Julia's optimized `sort!` directly. Otherwise ordering is an
insertion sort over adjacent transpositions, so every exchange is routed through
[`exchange_sign`](@ref) and the signs multiply out.
"""
function canonicalize_fields!(args::Vector{Field{S}}) where {S<:Statistics}
    if is_exchange_sign_free(S)
        sort!(args)
        return exchange_sign(S)
    end
    sign = Int8(1)
    @inbounds for i in 2:length(args)
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
