field_symbol(f::Field) = name(f)
field_symbol(f::Field{Boson}) = Symbol(string(name(f), is_classical(f) ? "ᶜ" : "ᴾ"))

function write_derivatives(io::IO, field::Field)
    for axis in derivative_multiindex(field)
        write(io, "∂", string(axis))
    end
    return nothing
end

function Base.show(io::IO, x::Field)
    reg = Int(regularisation(x))
    write_derivatives(io, x)
    write(io, string(field_symbol(x)))
    is_barred(x) && write(io, "̄")
    if reg == 1
        write(io, "⁺")
    elseif reg == -1
        write(io, "⁻")
    end
    return nothing
end

const show_brackets = Ref(true)
function Base.show(io::IO, x::QTerm)
    args = terms(x)
    isempty(args) && return nothing
    show_brackets[]::Bool && write(io, "(")
    show(io, args[1])
    f = SymbolicUtils.operation(x)
    for i in 2:length(args)
        show(io, f)
        show(io, args[i])
    end
    show_brackets[]::Bool && write(io, ")")
    return nothing
end

function Base.show(io::IO, x::QMul)
    if isempty(x.args_nc)
        show(io, x.arg_c)
        return nothing
    end
    if !SymbolicUtils._isone(x.arg_c)
        if SymbolicUtils._isone(-x.arg_c)
            write(io, "-")
        else
            print_number(io, x.arg_c)
            show(io, *)
        end
    end
    show_brackets[]::Bool && write(io, "(")
    show(io, x.args_nc[1])
    for i in 2:length(x.args_nc)
        show(io, *)
        show(io, x.args_nc[i])
    end
    show_brackets[]::Bool && write(io, ")")
    return nothing
end

function Base.show(io::IO, L::InteractionLagrangian)
    write(io, "Interaction Lagrangian with fields ")
    for (i, family) in enumerate(L.families)
        i > 1 && write(io, ", ")
        write(io, string(name(family)))
    end
    write(io, ":\n")
    lagrangian_terms = terms(L.lagrangian)
    if length(lagrangian_terms) == 1
        show(io, only(lagrangian_terms))
    else
        show(io, L.lagrangian)
    end
    return nothing
end

const T_LATEX = Union{QField,Diagrams,Diagram,Edge}
function Base.show(io::IO, ::MIME"text/latex", x::T_LATEX)
    write(io, "\$")
    write_latex(io, x)
    write(io, "\$")
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", L::InteractionLagrangian)
    write(io, "\$")
    write_latex(io, L)
    write(io, "\$")
    return nothing
end

const prop_type = Dict(
    PropagatorType.Retarded => "ᴿ",
    PropagatorType.Advanced => "ᴬ",
    PropagatorType.Keldysh => "ᴷ",
    PropagatorType.Spectral => "",
)
const reg_string = Dict(
    Regularisation.Plus => "⁺", Regularisation.Zero => "", Regularisation.Minus => "⁻"
)

function Base.show(io::IO, x::Edge)
    if !has_momenta(x)
        s = construct_position_basis(x)
    else
        s = construct_momentum_basis(x)
    end
    write(io, s)
    return nothing
end
function construct_position_basis(x::Edge)
    type = propagator_type(x)
    (out, in) = positions(x)
    (r2, r1) = regularisations(x)
    return string(
        is_spectral(type) ? "A" : "G",
        prop_type[type],
        "(",
        pos_string(out),
        reg_string[r2],
        ",",
        pos_string(in),
        reg_string[r1],
        ")",
    )
end
function construct_momentum_basis(x::Edge)
    type = propagator_type(x)
    m = repr(x.momenta)
    return string(is_spectral(type) ? "A" : "G", prop_type[type], "(", m, ")")
end

function Base.show(io::IO, d::Diagram)
    _contractions = contractions(d)
    l = length(_contractions)
    for idx in eachindex(_contractions)
        if idx == l
            show(io, _contractions[idx])
        else
            show(io, _contractions[idx])
            show(io, *)
        end
    end
    return nothing
end

_plain_negative_real(x::Real) = x < 0
function _plain_negative_real(x::Complex)
    return iszero(imag(x)) && real(x) < 0
end
_plain_negative_real(::Number) = false

function _show_key_value(io::IO, key, prefactor)
    print_number(io, prefactor)
    if !isempty(key)
        !SymbolicUtils._isone(prefactor) ? show(io, *) : write(io, "")
        show(io, key)
    end
    return nothing
end

function _show_dict_sum(io::IO, terms::AbstractDict)
    entries = collect(terms)
    sort!(entries; by=entry -> repr(first(entry)))
    for (i, (key, prefactor)) in enumerate(entries)
        if _plain_negative_real(prefactor)
            write(io, isone(i) ? "-" : " - ")
            _show_key_value(io, key, -prefactor)
        else
            i > 1 && write(io, " + ")
            _show_key_value(io, key, prefactor)
        end
    end
    return nothing
end

function Base.show(io::IO, ds::Diagrams)
    return _show_dict_sum(io, ds.diagrams)
end

function show_key(io::IO, terms::AbstractDict{K,C}, key::K) where {K,C}
    return _show_key_value(io, key, terms[key])
end

function print_number(io::IO, x::Real)
    if !SymbolicUtils._isone(x)
        show(io, x)
    end
    return nothing
end
function print_number(io::IO, x::Complex)
    if !SymbolicUtils._isone(x)
        if iszero(imag(x))
            show(io, real(x))
        else
            write(io, "(")
            show(io, x)
            write(io, ")")
        end
    end
    return nothing
end
function print_number(io::IO, x::Number)
    if !SymbolicUtils._isone(x)
        show(io, x)
    end
    return nothing
end

const _subscript_digits = Dict(
    '0' => "₀",
    '1' => "₁",
    '2' => "₂",
    '3' => "₃",
    '4' => "₄",
    '5' => "₅",
    '6' => "₆",
    '7' => "₇",
    '8' => "₈",
    '9' => "₉",
)

function _subscript_string(index::Integer)
    index >= 0 || throw(ArgumentError("display index must be non-negative"))
    return join(_subscript_digits[digit] for digit in string(index))
end

function pos_string(p::Position)
    if is_in(p)
        return "x₂"
    elseif is_out(p)
        return "x₁"
    else
        return "y" * _subscript_string(Int(p.index))
    end
end

function Base.show(io::IO, mime::MIME"text/plain", Σ::Union{SelfEnergy,DressedPropagator})
    Σ isa SelfEnergy ? print(io, "Self Energy:\n") : print(io, "Dressed Propagator:\n")
    write(io, "keldysh:  ")
    show(io, Σ.keldysh)
    write(io, "\nretarded: ")
    show(io, Σ.retarded)
    write(io, "\nadvanced: ")
    show(io, Σ.advanced)
    return nothing
end

function momentum_string(p::Momentum)
    if iszero(p.index)
        return "k"
    else
        return "q" * _subscript_string(Int(p.index))
    end
end

function Base.show(io::IO, p::Momentum)
    write(io, momentum_string(p))
    return nothing
end

function Base.show(io::IO, ms::Momenta)
    isempty(ms.prefactors) && return write(io, "0")
    wrote = false
    for (prefactor, momentum) in zip(ms.prefactors, ms.momenta)
        iszero(prefactor) && continue
        negative = prefactor < 0
        magnitude = abs(prefactor)
        if wrote
            write(io, negative ? " - " : " + ")
        elseif negative
            write(io, "-")
        end
        if magnitude != 1
            show(io, magnitude)
            write(io, "*")
        end
        show(io, momentum)
        wrote = true
    end
    wrote || write(io, "0")
    return nothing
end

function Base.show(io::IO, ds::BosonicDistributions)
    return _show_dict_sum(io, ds.terms)
end

function Base.show(io::IO, bds::BosonicDistributionTerm)
    for (i, ms) in enumerate(bds.momenta)
        i > 1 && write(io, "*")
        write(io, "F(")
        show(io, ms)
        write(io, ")")
    end
    return nothing
end

function Base.show(io::IO, ci::CollisionIntegral)
    write(io, "Collision integral:\n")
    show(io, ci.terms)
    return nothing
end
