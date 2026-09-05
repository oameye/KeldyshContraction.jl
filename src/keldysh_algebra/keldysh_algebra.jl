#########################
#      FieldFamily
#########################

"""Physical field identity shared by all Keldysh components."""
struct FieldFamily{S<:Statistics}
    name::Symbol
    indices::FieldIndices
end

FieldFamily{S}(name::Symbol) where {S<:Statistics} = FieldFamily{S}(name, FieldIndices())

name(f::FieldFamily) = f.name
statistics(::FieldFamily{S}) where {S<:Statistics} = S
field_indices(f::FieldFamily) = f.indices

function Base.isequal(a::FieldFamily{S}, b::FieldFamily{S}) where {S<:Statistics}
    return isequal(name(a), name(b)) && isequal(field_indices(a), field_indices(b))
end
Base.:(==)(a::FieldFamily{S}, b::FieldFamily{S}) where {S<:Statistics} = isequal(a, b)
function Base.isless(a::FieldFamily{S}, b::FieldFamily{S}) where {S<:Statistics}
    name(a) == name(b) || return isless(name(a), name(b))
    return isless(field_indices(a), field_indices(b))
end

#########################
#       Field{S}
#########################

"""Fundamental Keldysh path-integral field with statistics `S`."""
struct Field{S<:Statistics} <: QSym
    family::FieldFamily{S}
    orientation::Orientation.T
    keldysh::KeldyshIndex.T
    position::Position
    regularisation::Regularisation.T
end

"""Construct a Keldysh component of a field family."""
function Field(
    family::FieldFamily{S},
    keldysh::KeldyshIndex.T,
    orientation::Orientation.T=Orientation.Unbarred,
    reg::Regularisation.T=Regularisation.Zero,
    pos::Position=Bulk(),
) where {S<:Statistics}
    return Field{S}(family, orientation, keldysh, pos, reg)
end

Base.getindex(family::FieldFamily, keldysh::KeldyshIndex.T) = Field(family, keldysh)

FieldFamily(f::Field) = f.family
field_family(f::Field) = f.family
name(f::Field) = name(field_family(f))
field_indices(f::Field) = field_indices(field_family(f))

"""Copy a field, replacing the specified properties."""
function reconstruct(
    f::Field{S};
    family::FieldFamily{S}=field_family(f),
    orientation::Orientation.T=orientation(f),
    keldysh::KeldyshIndex.T=keldysh_index(f),
    position::Position=position(f),
    regularisation::Regularisation.T=regularisation(f),
) where {S<:Statistics}
    return Field(family, keldysh, orientation, regularisation, position)
end

statistics(::Field{S}) where {S<:Statistics} = S
orientation(f::Field) = f.orientation
keldysh_index(f::Field) = f.keldysh
regularisation(f::Field) = f.regularisation
position(f::Field) = f.position
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

"""Toggle the path-integral orientation of a field."""
function bar(f::Field)
    o = is_unbarred(f) ? Orientation.Barred : Orientation.Unbarred
    return reconstruct(f; orientation=o)
end

function Base.isequal(a::Field{S}, b::Field{S}) where {S}
    return isequal(field_family(a), field_family(b)) &&
           isequal(orientation(a), orientation(b)) &&
           isequal(keldysh_index(a), keldysh_index(b)) &&
           isequal(position(a), position(b)) &&
           isequal(regularisation(a), regularisation(b))
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

"""Exchange sign for statistics `S`."""
exchange_sign(::Type{Boson}) = Int8(1)

"""Exchange sign for two fields with statistics `S`."""
function exchange_sign(::Type{S}, ::Field{S}, ::Field{S}) where {S<:Statistics}
    return exchange_sign(S)
end

"""Return whether all field exchanges have sign `+1`."""
is_exchange_sign_free(::Type{S}) where {S<:Statistics} = isone(exchange_sign(S))

"""Sort fields in place and return the accumulated exchange sign."""
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

"""Declare physical field families."""
macro qfields(qs...)
    defs = map(qs) do q
        q isa Expr && q.head == :(::) && length(q.args) == 2 ||
            throw(ArgumentError("expected field declaration of the form name::Statistics"))

        fname, statistics_expr = q.args
        fname isa Symbol || throw(ArgumentError("field family name must be a symbol"))
        statistics_expr isa Symbol || throw(
            ArgumentError(
                "@qfields declares families only; construct components with family[Classical] or family[Quantum]"
            ),
        )

        return :($(esc(fname)) = FieldFamily{$(esc(statistics_expr))}($(QuoteNode(fname))))
    end

    names = map(q -> esc(q.args[1]), qs)
    return Expr(:block, defs..., :(tuple($(names...))))
end
