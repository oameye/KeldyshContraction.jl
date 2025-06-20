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
    (r2, r1) = regularisations(x)
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
function Base.show(io::IO, ds::Diagrams)
    diagrams = collect(keys(ds.diagrams))
    l = length(diagrams)
    for idx in eachindex(diagrams)
        if idx == l
            show_key(io, ds.diagrams, diagrams[idx])
        else
            show_key(io, ds.diagrams, diagrams[idx])
            write(io, " + ")
        end
    end
    return nothing
end
function show_key(io, terms::Dict, key)
    prefactor = terms[key]
    print_number(io, prefactor)
    if !isempty(key)
        !SymbolicUtils._isone(prefactor) ? show(io, *) : write(io, "")
        show(io, key)
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

function pos_string(p::Position)
    if is_in(p)
        return "x₂"
    elseif is_out(p)
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

function momentum_string(p::Momentum)
    if iszero(p.index)
        return "k"
    else
        return "q" * underscore_dict[p.index]
    end
end

function Base.show(io::IO, p::Momentum)
    write(io, momentum_string(p))
    return nothing
end

function Base.show(io::IO, ms::Momenta)
    if length(ms.prefactors) == 1 && iszero(ms.prefactors[1])
        return write(io, "0")
    end

    for (i, (p, m)) in enumerate(zip(ms.prefactors, ms.momenta))
        if i > 1
            if p < 0
                write(io, " - ")
            else
                write(io, " + ")
            end
        end
        show(io, m)
    end
    return nothing
end

function Base.show(io::IO, ds::BosonicDistributions)
    bd_terms = collect(keys(ds.terms))
    l = length(bd_terms)
    for idx in eachindex(bd_terms)
        if idx == l
            show_key(io, ds.terms, bd_terms[idx])
        else
            show_key(io, ds.terms, bd_terms[idx])
            write(io, " + ")
        end
    end
    return nothing
end

function Base.show(io::IO, bds::BosonicDistributionTerm)
    for (i, ms) in enumerate(bds.momenta)
        if i > 1
            write(io, "*")
        end
        write(io, string("F(", ms, ")"))
    end
    return nothing
end

function Base.show(io::IO, ci::CollisionIntegral)
    write(io, "Collision integral ")
    write(io, ":\n")
    show(io, ci.terms)
    return nothing
end
