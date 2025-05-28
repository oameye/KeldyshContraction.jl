# Base.show(io::IO, x::QSym) = write(io, name(x))
function Base.show(io::IO, x::Create)
    reg = Int(regularisation(x))
    if iszero(reg)
        s = string("̄", name(x))
    elseif reg == 1
        s = string("̄", name(x), "⁺")
    else
        s = string("̄", name(x), "⁻")
    end
    write(io, s)
    return nothing
end
function Base.show(io::IO, x::Destroy)
    reg = Int(regularisation(x))
    if iszero(reg)
        s = string(name(x))
    elseif reg == 1
        s = string(name(x), "⁺")
    else
        s = string(name(x), "⁻")
    end
    write(io, s)
    return nothing
end

const show_brackets = Ref(true)
function Base.show(io::IO, x::QTerm)
    show_brackets[]::Bool && write(io, "(")
    show(io, arguments(x)[1])
    f = SymbolicUtils.operation(x)
    for i in 2:length(arguments(x))
        show(io, f)
        show(io, arguments(x)[i])
    end
    show_brackets[]::Bool && write(io, ")")
    return nothing
end

function Base.show(io::IO, x::QMul)
    if !SymbolicUtils._isone(x.arg_c)
        print_number(io, x.arg_c)
        show(io, *)
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
    show(io, L.cfield)
    write(io, " and ")
    show(io, L.qfield)
    write(io, ":\n")
    show(io, L.lagrangian)
    return nothing
end

const T_LATEX = Union{<:QField,Diagrams,Diagram,Edge}
Base.show(io::IO, ::MIME"text/latex", x::T_LATEX) = write(io, latexify(x))
function Base.show(io::IO, ::MIME"text/latex", L::InteractionLagrangian)
    # write(io, "Interaction Lagrangian with fields ")
    # write(io, latexify(L.cfield))
    # write(io, " and ")
    # write(io, latexify(L.qfield))
    # println(io, ": \\newline")
    write(io, latexify(L.lagrangian))
    return nothing
end
# function Base.show(io::IO, ::MIME"text/latex", L::DressedPropagator)
#     return write(io,latexify([L.retarded,L.advanced, L.keldysh]))
# end
function Base.show(io::IO, x::Edge)
    prop_type = Dict(Retarded => "ᴿ", Advanced => "ᴬ", Keldysh => "ᴷ")

    reg_string = Dict(Plus => "⁺", Zero => "", Minus => "⁻")

    (out, in) = positions(x)
    (r2, r1) = regularisations(x)
    s = string(
        "G",
        prop_type[propagator_type(x)],
        "(",
        pos_string(out),
        reg_string[r2],
        ",",
        pos_string(in),
        reg_string[r1],
        ")",
    )
    write(io, s)
    return nothing
end
function Base.show(io::IO, d::Diagram)
    contractions = d.contractions
    l = length(contractions)
    for idx in eachindex(contractions)
        if idx == l
            show(io, contractions[idx])
        else
            show(io, contractions[idx])
            show(io, *)
        end
    end
    return nothing
end
function Base.show(io::IO, ds::Diagrams)
    diagrams = collect(keys(ds.diagrams))
    l = length(diagrams)
    for idx in eachindex(diagrams)
        if idx == l
            prefactor = ds.diagrams[diagrams[idx]]
            print_number(io, prefactor)
            !SymbolicUtils._isone(prefactor) ? show(io, *) : write(io, "")
            show(io, diagrams[idx])
        else
            prefactor = ds.diagrams[diagrams[idx]]
            print_number(io, prefactor)
            !SymbolicUtils._isone(prefactor) ? show(io, *) : write(io, "")
            show(io, diagrams[idx])
            write(io, " + ")
        end
    end
    return nothing
end

function print_number(io, x::Number)
    if !SymbolicUtils._isone(x)
        x = make_real(x)
        x isa Complex ? write(io, "(") : write(io, "")
        show(io, x)
        x isa Complex ? write(io, ")") : write(io, "")
    end
    return nothing
end

const underscore_dict = Dict(
    1 => "₁", 2 => "₂", 3 => "₃", 4 => "₄", 5 => "₅", 6 => "₆", 7 => "₇", 8 => "₈", 9 => "₉"
)

function pos_string(p)
    if p isa In
        return "x₂"
    elseif p isa Out
        return "x₁"
    else
        return "y" * underscore_dict[p.index]
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
