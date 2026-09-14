# Physics-facing displays for the staged kinetic compiler.
#
# These intentionally summarize the represented physics rather than exposing the storage
# layout of the intermediate IR. They are used by the REPL and by Documenter's executed
# examples; exact programmatic inspection remains available through semantic accessors.

function _show_stage_plain(io::IO, label::AbstractString, x)
    print(io, label, " (order=", order(x), ", statistics=", statistics(x), ", parameter=")
    show(io, parameters(x))
    print(io, ")")
    return nothing
end

function _show_component_counts(io::IO, x)
    print(
        io,
        "\n  components: K=",
        length(x.keldysh),
        ", R=",
        length(x.retarded),
        ", A=",
        length(x.advanced),
    )
    return nothing
end

# Extract the integer carried by Val without exposing it as user API.
_val_parameter(::Val{G}) where {G} = G

function _write_latex_parameter(io::IO, parameter::ParameterMonomial)
    if isempty(parameter.powers)
        write(io, "1")
        return nothing
    end
    for (i, power) in enumerate(parameter.powers)
        i > 1 && write(io, "\\,")
        write_latex_symbol(io, power.name)
        power.exponent == 1 || write(io, "^{", string(power.exponent), "}")
    end
    return nothing
end

function _write_latex_stage_metadata(io::IO, x)
    write(io, "\\qquad \\left[O=", string(order(x)), ",\\;\\mathcal P=")
    _write_latex_parameter(io, parameters(x))
    write(io, "\\right]")
    return nothing
end

function _latex_display(writer, io::IO)
    write(io, "\\[\n")
    writer()
    write(io, "\n\\]")
    return nothing
end

function Base.show(io::IO, ::MIME"text/latex", G::DressedPropagator)
    return _latex_display(io) do
        write(io, "G^{R,A,K}(x_1,x_2)")
        return _write_latex_stage_metadata(io, G)
    end
end

function Base.show(io::IO, ::MIME"text/plain", G::FourierDressedPropagator)
    _show_stage_plain(io, "Fourier dressed propagator", G)
    _show_component_counts(io, G)
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", G::FourierDressedPropagator)
    return _latex_display(io) do
        write(io, "G_F^{R,A,K}(p)")
        return _write_latex_stage_metadata(io, G)
    end
end

function Base.show(io::IO, ::MIME"text/plain", Σ::FourierSelfEnergy)
    _show_stage_plain(io, "Fourier 1PI self-energy", Σ)
    _show_component_counts(io, Σ)
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", Σ::FourierSelfEnergy)
    return _latex_display(io) do
        write(io, "\\Sigma_F^{R,A,K}(p)")
        return _write_latex_stage_metadata(io, Σ)
    end
end

function Base.show(io::IO, ::MIME"text/plain", G::WignerDressedPropagator)
    _show_stage_plain(io, "Wigner dressed propagator", G)
    print(io, "\n  gradient order: ", _val_parameter(gradient_order(G)))
    _show_component_counts(io, G)
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", G::WignerDressedPropagator)
    return _latex_display(io) do
        write(
            io,
            "G_W^{R,A,K}(X,k)\\big|_{\\nabla^{",
            string(_val_parameter(gradient_order(G))),
            "}}",
        )
        return _write_latex_stage_metadata(io, G)
    end
end

function Base.show(io::IO, ::MIME"text/plain", Σ::WignerSelfEnergy)
    _show_stage_plain(io, "Wigner 1PI self-energy", Σ)
    print(io, "\n  gradient order: ", _val_parameter(gradient_order(Σ)))
    _show_component_counts(io, Σ)
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", Σ::WignerSelfEnergy)
    return _latex_display(io) do
        write(
            io,
            "\\Sigma_W^{R,A,K}(X,k)\\big|_{\\nabla^{",
            string(_val_parameter(gradient_order(Σ))),
            "}}",
        )
        return _write_latex_stage_metadata(io, Σ)
    end
end

function Base.show(io::IO, ::MIME"text/plain", Σ::KineticSelfEnergy)
    _show_stage_plain(io, "Kinetic self-energy", Σ)
    _show_component_counts(io, Σ)
    print(io, "\n  representation: spectral/statistical, G^K = -i F A")
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", Σ::KineticSelfEnergy)
    return _latex_display(io) do
        write(io, "\\Sigma_{\\mathrm{kin}}^{R,A,K}[A,F],\\qquad G^K=-iFA")
        return _write_latex_stage_metadata(io, Σ)
    end
end

function Base.show(io::IO, ::MIME"text/plain", I::OffShellCollisionExpression)
    _show_stage_plain(io, "Off-shell Kadanoff-Baym collision", I)
    print(io, "\n  I_coll = i Σ^K - F_target(k) A_Σ,  A_Σ = i(Σ^R - Σ^A)")
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", I::OffShellCollisionExpression)
    return _latex_display(io) do
        write(
            io,
            "I_{\\mathrm{coll}}=i\\Sigma^K-F_{\\mathrm{target}}(k)A_{\\Sigma},\\qquad " *
            "A_{\\Sigma}=i(\\Sigma^R-\\Sigma^A)",
        )
        return _write_latex_stage_metadata(io, I)
    end
end

function Base.show(io::IO, ::MIME"text/plain", collision::SpectralDispersiveCollision)
    _show_stage_plain(io, "Spectral/dispersive collision", collision)
    print(io, "\n  G^R = D - iA/2,  G^A = D + iA/2")
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", collision::SpectralDispersiveCollision)
    return _latex_display(io) do
        write(
            io,
            "G^R=D-\\frac{i}{2}A,\\qquad G^A=D+\\frac{i}{2}A,\\qquad " *
            "I_{\\mathrm{coll}}=I_{\\mathrm{coll}}[D,A,F]",
        )
        return _write_latex_stage_metadata(io, collision)
    end
end

function Base.show(io::IO, ::MIME"text/plain", result::ReducedFrequencyCollision)
    _show_stage_plain(io, "Reduced causal-frequency collision", result)
    print(
        io,
        "\n  regular sectors: ",
        length(reduced_regular_terms(result)),
        "\n  causal blockers: ",
        length(reduced_blocked_terms(result)),
        "\n  unresolved Trotter states: ",
        length(reduced_trotter_terms(result)),
    )
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", result::ReducedFrequencyCollision)
    return _latex_display(io) do
        write(
            io,
            "I_{\\mathrm{coll}}\\longrightarrow " *
            "I_{\\mathrm{reg}}\\!\\left[\\delta(\\Delta E),\\operatorname{PV}\\!\\left(\\frac{1}{\\Delta E}\\right)\\right]" *
            "\\oplus I_{\\mathrm{blocked}}",
        )
        return write(
            io,
            "\\qquad (N_{\\mathrm{reg}},N_{\\mathrm{block}},N_{\\mathrm{T}})=(",
            string(length(reduced_regular_terms(result))),
            ",",
            string(length(reduced_blocked_terms(result))),
            ",",
            string(length(reduced_trotter_terms(result))),
            ")",
        )
    end
end

function Base.show(io::IO, ::MIME"text/plain", result::OccupationReducedExpression)
    _show_stage_plain(io, "Occupation-reduced collision", result)
    print(io, "\n  regular sectors: ", length(occupation_reduced_terms(result)))
    print(
        io,
        if statistics(result) === Boson
            "\n  F = 1 + 2n,  C_n = I/2"
        else
            "\n  F = 1 - 2n,  C_n = -I/2"
        end,
    )
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", result::OccupationReducedExpression)
    return _latex_display(io) do
        if statistics(result) === Boson
            write(io, "F=1+2n,\\qquad C_n=\\frac{1}{2}I_{\\mathrm{reg}}")
        else
            write(io, "F=1-2n,\\qquad C_n=-\\frac{1}{2}I_{\\mathrm{reg}}")
        end
        return write(io, "\\qquad N_{\\mathrm{sec}}=", string(length(result)))
    end
end

function Base.show(io::IO, ::MIME"text/plain", result::LoopQuotientedExpression)
    _show_stage_plain(io, "Loop-quotiented collision", result)
    print(io, "\n  canonical sectors: ", length(result))
    print(io, "\n  dummy loops quotiented by signed permutations")
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", result::LoopQuotientedExpression)
    return _latex_display(io) do
        return write(
            io,
            "\\{q_i\\}\\,/\\,\\big(q_i\\sim \\pm q_{\\pi(i)}\\big)" *
            "\\qquad\\Longrightarrow\\qquad N_{\\mathrm{sec}}=",
            string(length(result)),
        )
    end
end

function Base.show(io::IO, ::MIME"text/plain", kernel::CollisionKernel)
    _show_stage_plain(io, "Collision kernel", kernel)
    print(io, "\n  canonical sectors: ", length(kernel))
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", kernel::CollisionKernel)
    return _latex_display(io) do
        write(
            io,
            "C_n(k)=\\sum_{\\alpha=1}^{",
            string(length(kernel)),
            "}\\mathcal K_{\\alpha}(k,\\{q\\})\\," *
            "P_{\\alpha}[n]\\,\\mathcal S_{\\alpha}",
        )
        return _write_latex_stage_metadata(io, kernel)
    end
end
