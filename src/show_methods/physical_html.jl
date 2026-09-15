# Renderer-safe HTML displays for the physical bosonic and fermionic kinetic pipeline.
#
# Documenter prefers `text/html` over `text/latex`. Splitting large symbolic sums into
# separate math containers avoids feeding one oversized/alignment-heavy TeX expression
# to the browser math engine while preserving the package's `text/plain` and `text/latex`
# representations.

const _PHYSICAL_LATEX_SUM_BREAK = "\\\\\n&{}"

function _html_escape_math(text::AbstractString)
    return replace(text, '&' => "&amp;", '<' => "&lt;", '>' => "&gt;")
end

function _write_html_math(io::IO, latex::AbstractString)
    write(io, "<div class=\"math-container\">\\[", _html_escape_math(latex), "\\]</div>")
    return nothing
end

function _physical_latex_terms(expression::AbstractString)
    return split(expression, _PHYSICAL_LATEX_SUM_BREAK; keepempty=false)
end

function _write_html_sum(io::IO, lhs::AbstractString, expression::AbstractString)
    terms = _physical_latex_terms(expression)
    if isempty(terms)
        return _write_html_math(io, lhs * "=0")
    end

    _write_html_math(io, lhs * "={}" * first(terms))
    continuation = "\\phantom{" * lhs * "={}}"
    for term in Iterators.drop(terms, 1)
        _write_html_math(io, continuation * term)
    end
    return nothing
end

function _physical_component_latex(component, parameter::ParameterMonomial)
    if component isa Diagrams
        return _coordinate_component_string(component, parameter, true)
    elseif component isa Union{FourierDiagrams,WignerDiagrams}
        return _routed_component_string(component, parameter, true)
    end
    return _kinetic_component_string(component, parameter, true)
end

function _html_three_components(
    io::IO, object, symbols::NTuple{3,String}, argument::AbstractString
)
    parameter = parameters(object)
    components = (object.keldysh, object.retarded, object.advanced)
    write(io, "<div class=\"kc-physics-display\">")
    for (symbol, component) in zip(symbols, components)
        lhs = symbol * "\\!\\left(" * argument * "\\right)"
        _write_html_sum(io, lhs, _physical_component_latex(component, parameter))
    end
    write(io, "</div>")
    return nothing
end

function Base.show(
    io::IO, ::MIME"text/html", G::DressedPropagator{C,S,O,E1,E2}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2}
    return _html_three_components(io, G, ("G^K", "G^R", "G^A"), "x_1,x_2")
end

function Base.show(
    io::IO, ::MIME"text/html", Σ::SelfEnergy{C,S,O,E1,E2}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2}
    return _html_three_components(io, Σ, ("\\Sigma^K", "\\Sigma^R", "\\Sigma^A"), "x_1,x_2")
end

function Base.show(
    io::IO, ::MIME"text/html", G::FourierDressedPropagator{C,S,O,E1,E2}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2}
    return _html_three_components(io, G, ("G_F^K", "G_F^R", "G_F^A"), "k")
end

function Base.show(
    io::IO, ::MIME"text/html", Σ::FourierSelfEnergy{C,S,O,E1,E2}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2}
    return _html_three_components(io, Σ, ("\\Sigma_F^K", "\\Sigma_F^R", "\\Sigma_F^A"), "k")
end

function Base.show(
    io::IO, ::MIME"text/html", G::WignerDressedPropagator{C,S,O,E1,E2,GOrder,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,GOrder,Ctx<:AbstractWignerContext}
    return _html_three_components(io, G, ("G_W^K", "G_W^R", "G_W^A"), "X,k")
end

function Base.show(
    io::IO, ::MIME"text/html", Σ::WignerSelfEnergy{C,S,O,E1,E2,GOrder,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,GOrder,Ctx<:AbstractWignerContext}
    return _html_three_components(
        io, Σ, ("\\Sigma_W^K", "\\Sigma_W^R", "\\Sigma_W^A"), "X,k"
    )
end

function Base.show(
    io::IO, ::MIME"text/html", Σ::KineticSelfEnergy{C,S,O,E1,E2,GOrder,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,GOrder,Ctx<:AbstractWignerContext}
    return _html_three_components(
        io,
        Σ,
        (
            "\\Sigma_{\\mathrm{kin}}^K",
            "\\Sigma_{\\mathrm{kin}}^R",
            "\\Sigma_{\\mathrm{kin}}^A",
        ),
        "k",
    )
end

function Base.show(
    io::IO,
    ::MIME"text/html",
    collision::OffShellCollisionExpression{C,S,O,E1,E2,GOrder,Ctx},
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,GOrder,Ctx<:AbstractWignerContext}
    write(io, "<div class=\"kc-physics-display\">")
    _write_html_sum(io, "I_{\\mathrm{coll}}(k)", _off_shell_string(collision, true))
    write(io, "</div>")
    return nothing
end

function Base.show(
    io::IO,
    ::MIME"text/html",
    collision::SpectralDispersiveCollision{C,S,O,E1,E2,GOrder,Ctx},
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,GOrder,Ctx<:AbstractWignerContext}
    write(io, "<div class=\"kc-physics-display\">")
    _write_html_sum(
        io, "I_{\\mathrm{coll}}(k)", _spectral_dispersive_string(collision, true)
    )
    write(io, "</div>")
    return nothing
end

function Base.show(
    io::IO, ::MIME"text/html", result::ReducedFrequencyCollision{C,S,O,E1,E2,GOrder,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,GOrder,Ctx<:AbstractWignerContext}
    write(io, "<div class=\"kc-physics-display\">")
    _write_html_sum(
        io,
        "I_{\\mathrm{reg}}(k)",
        _physical_collision_string(reduced_regular_terms(result), true),
    )
    _write_html_math(io, _blocked_summary(result, true))
    write(io, "</div>")
    return nothing
end

function Base.show(
    io::IO, ::MIME"text/html", result::OccupationReducedExpression{C,S,O,GOrder,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,GOrder,Ctx<:AbstractWignerContext}
    write(io, "<div class=\"kc-physics-display\">")
    _write_html_sum(
        io,
        "C_n^{\\mathrm{reg}}(k)",
        _physical_collision_string(occupation_reduced_terms(result), true),
    )
    write(io, "</div>")
    return nothing
end

function Base.show(
    io::IO, ::MIME"text/html", result::LoopQuotientedExpression{C,S,O,GOrder,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,GOrder,Ctx<:AbstractWignerContext}
    write(io, "<div class=\"kc-physics-display\">")
    _write_html_sum(
        io,
        "C_n^{\\mathrm{quot}}(k)",
        _physical_collision_string(loop_quotient_terms(result), true),
    )
    write(io, "</div>")
    return nothing
end

function Base.show(
    io::IO, ::MIME"text/html", kernel::CollisionKernel{C,S,O,GOrder,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,GOrder,Ctx<:AbstractWignerContext}
    write(io, "<div class=\"kc-physics-display\">")
    _write_html_sum(
        io, "C_n(k)", _physical_collision_string(collision_kernel_terms(kernel), true)
    )
    write(io, "</div>")
    return nothing
end
