#########################
# Destroy and Create
#########################

"""
    Destroy <: QSym

Bosonic field representing the quantum field annihilation operator.
"""
struct Destroy{P} <: QSym
    name::Symbol
    contour::KeldyshContour
    position::P
    regularisation::Regularisation
    function Destroy(
        name::Symbol,
        contour::KeldyshContour,
        reg::Regularisation=Zero,
        pos::AbstractPosition=Bulk(),
    )
        return new{typeof(pos)}(name, contour, pos, reg)
    end
end
# TODO: now we can dispatch on Bulk and In as they are different structs.
# Instead, we want to dispatch on the position of the field.
# Because, otherwise later you get in problem when trying to make Diagram concrete

"""
    Create <: QSym

Bosonic field representing the quantum field creation operator.
"""
struct Create{P} <: QSym
    name::Symbol
    contour::KeldyshContour
    position::P
    regularisation::Regularisation
    function Create(
        name::Symbol,
        contour::KeldyshContour,
        reg::Regularisation=Zero,
        pos::AbstractPosition=Bulk();
    )
        return new{typeof(pos)}(name, contour, pos, reg)
    end
end

for T in (:Create, :Destroy)
    @eval Base.isequal(a::$T, b::$T) =
        isequal(name(a), name(b)) &&
        isequal(position(a), position(b)) &&
        isequal(contour(a), contour(b)) &&
        isequal(regularisation(a), regularisation(b))
end

for f in [:Destroy, :Create]
    @eval function (ff::$f)(pos::AbstractPosition)
        return $(f)(name(ff), contour(ff), regularisation(ff), pos)
    end
    @eval function (ff::$f)(reg::Regularisation)
        return $(f)(name(ff), contour(ff), reg, position(ff))
    end

    @eval regularisation(ϕ::$(f)) = ϕ.regularisation
    @eval contour(ϕ::$(f)) = ϕ.contour
    @eval position(ϕ::$(f)) = ϕ.position
    @eval isbulk(ϕ::$(f)) = isbulk(position(ϕ))
    @eval is_in(ϕ::$(f)) = is_in(position(ϕ))
    @eval is_out(ϕ::$(f)) = is_out(position(ϕ))

    @eval set_reg_to_zero(ϕ::$(f)) = $(f)(name(ϕ), contour(ϕ), Zero, position(ϕ))
end
function set_reg_to_zero!(v::Vector{QSym})
    for (i, f) in pairs(v)
        v[i] = set_reg_to_zero(f)
    end
end
function contour_integers(v::Vector{QSym})
    Int[Int(contour(x)) for x in v]
end
ladder(::Destroy) = 0
ladder(::Create) = 1

"""
    adjoint(op::Destroy)

Adjoint of the operator [`Destroy`](@ref) annihilation field constructing the corresponding creation field [`Create`](@ref).
"""
function Base.adjoint(op::Destroy)
    return Create(name(op), contour(op), regularisation(op), position(op))
end

"""
    adjoint(op::Create)

Adjoint of the [`Create`](@ref) creation field constructing the corresponding annihilation field [`Destroy`](@ref).
"""
function Base.adjoint(op::Create)
    return Destroy(name(op), contour(op), regularisation(op), position(op))
end

# reverse normal ordered
"Defines the reverse normal order for the creation and annihilation operators."
function Base.:*(a::Create, b::Destroy)
    pos_a = position(a)
    pos_b = position(b)
    if pos_a < pos_b
        return QMul(1, [a, b])
    else
        return QMul(1, [b, a])
    end
end

"""
    is_creation(x::QSym)

Returns if the field `x` is a creation operator.
"""
is_creation(x::Create) = true
is_creation(x::Destroy) = false

"""
    is_annihilation(x::QSym)

Returns if the field `x` is an annihilation operator.
"""
is_annihilation(x::Destroy) = true
is_annihilation(x::Create) = false

"""
    @qfields

Convenience macro for the construction of fields.

Examples
========
```
julia> using KeldyshContraction: Classical, Quantum

julia> @qnumbers ψ::Destroy(Classical)
(ψ,)
```
"""
macro qfields(qs...)
    defs = map(qs) do q
        nf = _name_field(q)
        name, field_expr = nf.name, nf.field_expr

        # Extract field type (Create or Destroy) and args
        field_type = field_expr.args[1]
        field_args = field_expr.args[2:end]

        # Build the field construction expression
        field_construction = Expr(
            :call,
            esc(field_type),
            Expr(:quote, name),  # Pass the name as a quoted symbol
            map(esc, field_args)..., # Escape all other args
        )

        # Define the field variable
        :($(esc(name)) = $(field_construction))
    end

    # Return both the definitions and the tuple of names
    return Expr(:block, defs..., :(tuple($(map(q -> esc(_name_field(q).name), qs)...))))
end

# Helper function to extract name and field constructor from expression
function _name_field(expr)
    @assert expr isa Expr && expr.head == :(::) "Expected expression of form name::Field(args...)"
    name = expr.args[1]
    @assert name isa Symbol "Left side of :: must be a symbol"

    field_expr = expr.args[2]
    @assert field_expr isa Expr && field_expr.head == :call "Right side of :: must be a field constructor call"

    return (name=name, field_expr=field_expr)
end
