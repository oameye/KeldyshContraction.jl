const LATEX_SYMBOLS = Dict{Symbol,String}(
    :α => "\\alpha",
    :β => "\\beta",
    :γ => "\\gamma",
    :δ => "\\delta",
    :ϵ => "\\epsilon",
    :ε => "\\varepsilon",
    :ζ => "\\zeta",
    :η => "\\eta",
    :θ => "\\theta",
    :ϑ => "\\vartheta",
    :ι => "\\iota",
    :κ => "\\kappa",
    :λ => "\\lambda",
    :μ => "\\mu",
    :ν => "\\nu",
    :ξ => "\\xi",
    :π => "\\pi",
    :ϖ => "\\varpi",
    :ρ => "\\rho",
    :ϱ => "\\varrho",
    :σ => "\\sigma",
    :ς => "\\varsigma",
    :τ => "\\tau",
    :υ => "\\upsilon",
    :ϕ => "\\phi",
    :φ => "\\varphi",
    :χ => "\\chi",
    :ψ => "\\psi",
    :ω => "\\omega",
    :Γ => "\\Gamma",
    :Δ => "\\Delta",
    :Θ => "\\Theta",
    :Λ => "\\Lambda",
    :Ξ => "\\Xi",
    :Π => "\\Pi",
    :Σ => "\\Sigma",
    :Υ => "\\Upsilon",
    :Φ => "\\Phi",
    :Ψ => "\\Psi",
    :Ω => "\\Omega",
)

function write_latex_symbol(io::IO, symbol::Symbol)
    if haskey(LATEX_SYMBOLS, symbol)
        write(io, LATEX_SYMBOLS[symbol])
    else
        write(io, string(symbol))
    end
    return nothing
end

function write_latex_regularisation(io::IO, field::Field)
    reg = Int(regularisation(field))
    if reg == 1
        write(io, "^+")
    elseif reg == -1
        write(io, "^{-}")
    end
    return nothing
end

function write_latex_field(io::IO, field::Field{Boson}, ::Val{:standalone})
    if is_barred(field)
        write(io, "\\bar{", string(name(field)), is_classical(field) ? "ᶜ" : "ᴾ", "}")
    else
        write_latex_symbol(io, name(field))
        write(io, is_classical(field) ? "^c" : "^P")
    end
    write_latex_regularisation(io, field)
    return nothing
end

function write_latex_field(io::IO, field::Field{Boson}, ::Val{:expression})
    if is_barred(field)
        write(io, "\\bar{")
        write_latex_symbol(io, name(field))
        write(io, is_classical(field) ? "^c}" : "^P}")
    else
        write_latex_symbol(io, name(field))
        write(io, is_classical(field) ? "^c" : "^P")
    end
    write_latex_regularisation(io, field)
    return nothing
end

function write_latex(io::IO, field::Field{Boson})
    write_latex_field(io, field, Val(:standalone))
    return nothing
end

function write_latex_number(io::IO, x::Rational)
    if denominator(x) == 1
        show(io, numerator(x))
    else
        write(io, "\\frac{", string(numerator(x)), "}{", string(denominator(x)), "}")
    end
    return nothing
end
function write_latex_number(io::IO, x::Real)
    show(io, x)
    return nothing
end
function write_latex_number(io::IO, x::Complex)
    re = real(x)
    im_part = imag(x)
    if iszero(im_part)
        write_latex_number(io, re)
    elseif iszero(re)
        if isone(im_part)
            write(io, "i")
        elseif isone(-im_part)
            write(io, "-i")
        else
            write_latex_number(io, im_part)
            write(io, " i")
        end
    else
        write(io, "(")
        write_latex_number(io, re)
        if im_part < 0
            write(io, " - ")
            write_latex_number(io, -im_part)
        else
            write(io, " + ")
            write_latex_number(io, im_part)
        end
        write(io, " i)")
    end
    return nothing
end
function write_latex_number(io::IO, x::Number)
    show(io, x)
    return nothing
end

function write_latex(io::IO, term::QMul)
    coefficient = term.arg_c
    fields = term.args_nc
    if isempty(fields)
        write_latex_number(io, coefficient)
        return nothing
    end

    if !SymbolicUtils._isone(coefficient)
        if SymbolicUtils._isone(-coefficient)
            write(io, "- ")
        else
            write_latex_number(io, coefficient)
            write(io, " ")
        end
    end
    for (i, field) in enumerate(fields)
        i > 1 && write(io, " ")
        write_latex_field(io, field, Val(:expression))
    end
    return nothing
end

function write_latex(io::IO, sum::QAdd)
    for (i, term) in enumerate(sum.arguments)
        i > 1 && write(io, " + ")
        write_latex(io, term)
    end
    return nothing
end

function write_latex_position(io::IO, position::Position, reg::Regularisation.T)
    if is_in(position)
        write(io, "x_2")
    elseif is_out(position)
        write(io, "x_1")
    else
        write(io, "y_", string(Int(index(position))))
    end
    if reg === Regularisation.Plus
        write(io, "^+")
    elseif reg === Regularisation.Minus
        write(io, "^-")
    end
    return nothing
end

function write_latex(io::IO, momentum::Momentum)
    if iszero(momentum.index)
        write(io, "k")
    else
        write(io, "q_", string(Int(momentum.index)))
    end
    return nothing
end

function write_latex(io::IO, momenta::Momenta)
    if length(momenta.prefactors) == 1 && iszero(momenta.prefactors[1])
        write(io, "0")
        return nothing
    end
    for (i, (prefactor, momentum)) in enumerate(zip(momenta.prefactors, momenta.momenta))
        if i > 1
            write(io, prefactor < 0 ? " - " : " + ")
        elseif prefactor < 0
            write(io, "-")
        end
        write_latex(io, momentum)
    end
    return nothing
end

function write_latex(io::IO, edge::Edge)
    type = propagator_type(edge)
    write(io, is_spectral(type) ? "A" : "G")
    if is_retarded(type)
        write(io, "^R")
    elseif is_advanced(type)
        write(io, "^A")
    elseif is_keldysh(type)
        write(io, "^K")
    end
    write(io, "\\left( ")
    if has_momenta(edge)
        write_latex(io, edge.momenta)
    else
        out, in = positions(edge)
        r2, r1 = regularisations(edge)
        write_latex_position(io, out, r2)
        write(io, ", ")
        write_latex_position(io, in, r1)
    end
    write(io, " \\right)")
    return nothing
end

function write_latex(io::IO, diagram::Diagram)
    for (i, edge) in enumerate(contractions(diagram))
        i > 1 && write(io, " ")
        write_latex(io, edge)
    end
    return nothing
end

function write_latex(io::IO, diagrams::Diagrams)
    entries = collect(diagrams.diagrams)
    for (i, (diagram, coefficient)) in enumerate(entries)
        i > 1 && write(io, " + ")
        if !SymbolicUtils._isone(coefficient)
            write_latex_number(io, coefficient)
            isempty(diagram) || write(io, " ")
        end
        write_latex(io, diagram)
    end
    return nothing
end

function write_latex(io::IO, L::InteractionLagrangian)
    write_latex(io, L.lagrangian)
    return nothing
end
