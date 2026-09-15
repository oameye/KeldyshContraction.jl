# Complete physics-facing displays for the ordinary bosonic and fermionic compiler paths.
#
# The generic stage summaries in `kinetic_pipeline.jl` remain a fallback for future statistics
# extensions. The physical statistics supported by the package render the represented expression
# itself at every stage of the derivation.

const _PhysicalDisplayStatistics = Union{Boson,Fermion}

function _render_with(writer, value)
    io = IOBuffer()
    writer(io, value)
    return String(take!(io))
end

_plain_render(value) = sprint(show, value)
_latex_render(value) = _render_with(write_latex, value)

function _coordinate_component_string(
    diagrams::Diagrams, parameter::ParameterMonomial, latex::Bool
)
    rendered = String[]
    parameter_text = _parameter_string(parameter, latex)
    for (diagram, coefficient) in diagrams
        diagram_text = latex ? _latex_render(diagram) : _plain_render(diagram)
        push!(
            rendered, _term_string(coefficient, String[parameter_text, diagram_text], latex)
        )
    end
    sort!(rendered)
    return _sum_string(rendered, latex)
end

function _edge_shift(edge::Edge)
    out_reg, in_reg = regularisations(edge)
    return Int(out_reg) - Int(in_reg)
end

function _propagator_display_parts(type::PropagatorType.T)
    is_spectral(type) && return ("A", nothing)
    is_retarded(type) && return ("G", "R")
    is_advanced(type) && return ("G", "A")
    is_keldysh(type) && return ("G", "K")
    return ("G", nothing)
end

function _line_display_string(
    symbol::String,
    component::Union{Nothing,String},
    family_text::String,
    momentum_text::String,
    latex::Bool,
)
    if latex
        decorated = symbol * "_{" * family_text * "}"
        component === nothing || (decorated *= "^{" * component * "}")
        return decorated * "\\!\\left(" * momentum_text * "\\right)"
    end
    decorated = component === nothing ? symbol : symbol * "^" * component
    return decorated * "_" * family_text * "(" * momentum_text * ")"
end

function _routed_line_string(
    edge::Edge,
    routed::LinearMomentum,
    basis::MomentumBasis,
    external::MomentumVariable,
    latex::Bool,
)
    family = field_family(first(fields(edge)))
    type = propagator_type(edge)
    family_text = _regulated_family_string(family, _edge_shift(edge), latex)
    momentum_text = _linear_momentum_string(routed, basis, external, latex)
    symbol, component = _propagator_display_parts(type)
    return _line_display_string(symbol, component, family_text, momentum_text, latex)
end

_display_external(graph::FourierDiagram) = momentum_basis(graph)[1]
_display_external(graph::WignerDiagram) = external_wigner_momentum(graph)

function _routed_graph_terms(
    graph, contributions, parameter::ParameterMonomial, latex::Bool
)
    rendered = String[]
    basis = momentum_basis(graph)
    external = _display_external(graph)
    parameter_text = _parameter_string(parameter, latex)
    measure_text = _measure_string(basis, external, latex)
    coordinate = coordinate_diagram(graph)
    line_factors = String[
        _routed_line_string(edge, routed, basis, external, latex) for
        (edge, routed) in zip(contractions(coordinate), edge_momenta(graph))
    ]

    for contribution in contributions
        for (monomial, kinematic_coefficient) in contribution.kinematic
            coefficient = contribution.coefficient * kinematic_coefficient
            iszero(coefficient) && continue
            kinematic_text = _momentum_monomial_string(monomial, basis, external, latex)
            factors = vcat(
                String[parameter_text, measure_text, kinematic_text], line_factors
            )
            push!(rendered, _term_string(coefficient, factors, latex))
        end
    end
    return rendered
end

function _routed_component_string(diagrams, parameter::ParameterMonomial, latex::Bool)
    rendered = String[]
    for (graph, contributions) in diagrams
        append!(rendered, _routed_graph_terms(graph, contributions, parameter, latex))
    end
    sort!(rendered)
    return _sum_string(rendered, latex)
end

function _kinetic_line_string(
    line::KineticLine, basis::MomentumBasis, external::MomentumVariable, latex::Bool
)
    kind = kinetic_line_kind(line)
    symbol, component = if kind === KineticSpectral
        ("A", nothing)
    elseif kind === KineticRetarded
        ("G", "R")
    else
        ("G", "A")
    end
    family_text = _regulated_family_string(
        line.family, Int(regularisation_shift(line)), latex
    )
    momentum_text = _linear_momentum_string(momentum(line), basis, external, latex)
    line_text = _line_display_string(symbol, component, family_text, momentum_text, latex)

    statistical_weight(line) === DistributionWeight || return line_text
    distribution = _distribution_atom_string(
        "F", line.family, momentum(line), basis, external, latex
    )
    return latex ? distribution * "\\," * line_text : distribution * " " * line_text
end

function _kinetic_expression_terms(
    expression::KineticExpression,
    parameter::ParameterMonomial,
    latex::Bool;
    external_distribution::Union{Nothing,FieldFamily}=nothing,
)
    rendered = String[]
    parameter_text = _parameter_string(parameter, latex)
    for (term, source_coefficient) in expression
        basis = momentum_basis(term)
        external = external_wigner_momentum(term)
        measure_text = _measure_string(basis, external, latex)
        line_factors = String[
            _kinetic_line_string(line, basis, external, latex) for
            line in kinetic_lines(term)
        ]
        if external_distribution !== nothing
            external_slot = _external_basis_index(basis, external)
            external_momentum = basis_momentum(basis, external_slot)
            pushfirst!(
                line_factors,
                _distribution_atom_string(
                    "F", external_distribution, external_momentum, basis, external, latex
                ),
            )
        end
        for (monomial, kinematic_coefficient) in kinematic_factor(term)
            coefficient = source_coefficient * kinematic_coefficient
            iszero(coefficient) && continue
            kinematic_text = _momentum_monomial_string(monomial, basis, external, latex)
            factors = vcat(
                String[parameter_text, measure_text, kinematic_text], line_factors
            )
            push!(rendered, _term_string(coefficient, factors, latex))
        end
    end
    sort!(rendered)
    return rendered
end

function _kinetic_component_string(
    expression::KineticExpression, parameter::ParameterMonomial, latex::Bool
)
    return _sum_string(_kinetic_expression_terms(expression, parameter, latex), latex)
end

function _off_shell_string(collision::OffShellCollisionExpression, latex::Bool)
    rendered = _kinetic_expression_terms(collision.offset, parameters(collision), latex)
    append!(
        rendered,
        _kinetic_expression_terms(
            collision.distribution_coefficient,
            parameters(collision),
            latex;
            external_distribution=target_family(collision),
        ),
    )
    sort!(rendered)
    return _sum_string(rendered, latex)
end

function _plain_three_components(
    io::IO, title::AbstractString, object, symbols::NTuple{3,String}
)
    _show_stage_plain(io, title, object)
    parameter = parameters(object)
    components = (object.keldysh, object.retarded, object.advanced)
    for (symbol, component) in zip(symbols, components)
        write(io, "\n  ", symbol, " = ")
        if component isa Diagrams
            write(io, _coordinate_component_string(component, parameter, false))
        elseif component isa Union{FourierDiagrams,WignerDiagrams}
            write(io, _routed_component_string(component, parameter, false))
        else
            write(io, _kinetic_component_string(component, parameter, false))
        end
    end
    return nothing
end

function _latex_three_components(
    io::IO, object, symbols::NTuple{3,String}, argument::AbstractString
)
    parameter = parameters(object)
    components = (object.keldysh, object.retarded, object.advanced)
    return _latex_display(io) do
        write(io, "\\begin{aligned}\n")
        for i in eachindex(symbols)
            i > 1 && write(io, "\\\\\n")
            write(io, symbols[i], "\\!\\left(", argument, "\\right)&={}")
            component = components[i]
            if component isa Diagrams
                write(io, _coordinate_component_string(component, parameter, true))
            elseif component isa Union{FourierDiagrams,WignerDiagrams}
                write(io, _routed_component_string(component, parameter, true))
            else
                write(io, _kinetic_component_string(component, parameter, true))
            end
        end
        write(io, "\n\\end{aligned}")
        return nothing
    end
end

function Base.show(
    io::IO, ::MIME"text/plain", G::DressedPropagator{C,S,O,E1,E2}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2}
    return _plain_three_components(io, "Dressed propagator", G, ("G^K", "G^R", "G^A"))
end

function Base.show(
    io::IO, ::MIME"text/latex", G::DressedPropagator{C,S,O,E1,E2}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2}
    return _latex_three_components(io, G, ("G^K", "G^R", "G^A"), "x_1,x_2")
end

function Base.show(
    io::IO, ::MIME"text/plain", Σ::SelfEnergy{C,S,O,E1,E2}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2}
    return _plain_three_components(io, "1PI self-energy", Σ, ("Σ^K", "Σ^R", "Σ^A"))
end

function Base.show(
    io::IO, ::MIME"text/latex", Σ::SelfEnergy{C,S,O,E1,E2}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2}
    return _latex_three_components(
        io, Σ, ("\\Sigma^K", "\\Sigma^R", "\\Sigma^A"), "x_1,x_2"
    )
end

function Base.show(
    io::IO, ::MIME"text/plain", G::FourierDressedPropagator{C,S,O,E1,E2}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2}
    return _plain_three_components(
        io, "Fourier dressed propagator", G, ("G_F^K", "G_F^R", "G_F^A")
    )
end

function Base.show(
    io::IO, ::MIME"text/latex", G::FourierDressedPropagator{C,S,O,E1,E2}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2}
    return _latex_three_components(io, G, ("G_F^K", "G_F^R", "G_F^A"), "k")
end

function Base.show(
    io::IO, ::MIME"text/plain", Σ::FourierSelfEnergy{C,S,O,E1,E2}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2}
    return _plain_three_components(
        io, "Fourier 1PI self-energy", Σ, ("Σ_F^K", "Σ_F^R", "Σ_F^A")
    )
end

function Base.show(
    io::IO, ::MIME"text/latex", Σ::FourierSelfEnergy{C,S,O,E1,E2}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2}
    return _latex_three_components(
        io, Σ, ("\\Sigma_F^K", "\\Sigma_F^R", "\\Sigma_F^A"), "k"
    )
end

function Base.show(
    io::IO, ::MIME"text/plain", G::WignerDressedPropagator{C,S,O,E1,E2,GOrder,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,GOrder,Ctx<:AbstractWignerContext}
    return _plain_three_components(
        io, "Wigner dressed propagator", G, ("G_W^K", "G_W^R", "G_W^A")
    )
end

function Base.show(
    io::IO, ::MIME"text/latex", G::WignerDressedPropagator{C,S,O,E1,E2,GOrder,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,GOrder,Ctx<:AbstractWignerContext}
    return _latex_three_components(io, G, ("G_W^K", "G_W^R", "G_W^A"), "X,k")
end

function Base.show(
    io::IO, ::MIME"text/plain", Σ::WignerSelfEnergy{C,S,O,E1,E2,GOrder,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,GOrder,Ctx<:AbstractWignerContext}
    return _plain_three_components(
        io, "Wigner 1PI self-energy", Σ, ("Σ_W^K", "Σ_W^R", "Σ_W^A")
    )
end

function Base.show(
    io::IO, ::MIME"text/latex", Σ::WignerSelfEnergy{C,S,O,E1,E2,GOrder,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,GOrder,Ctx<:AbstractWignerContext}
    return _latex_three_components(
        io, Σ, ("\\Sigma_W^K", "\\Sigma_W^R", "\\Sigma_W^A"), "X,k"
    )
end

function Base.show(
    io::IO, ::MIME"text/plain", Σ::KineticSelfEnergy{C,S,O,E1,E2,GOrder,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,GOrder,Ctx<:AbstractWignerContext}
    return _plain_three_components(
        io, "Kinetic self-energy", Σ, ("Σ_kin^K", "Σ_kin^R", "Σ_kin^A")
    )
end

function Base.show(
    io::IO, ::MIME"text/latex", Σ::KineticSelfEnergy{C,S,O,E1,E2,GOrder,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,GOrder,Ctx<:AbstractWignerContext}
    return _latex_three_components(
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
    ::MIME"text/plain",
    collision::OffShellCollisionExpression{C,S,O,E1,E2,GOrder,Ctx},
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,GOrder,Ctx<:AbstractWignerContext}
    _show_stage_plain(io, "Off-shell Kadanoff-Baym collision", collision)
    write(io, "\n  I_coll(k) = ", _off_shell_string(collision, false))
    return nothing
end

function Base.show(
    io::IO,
    ::MIME"text/latex",
    collision::OffShellCollisionExpression{C,S,O,E1,E2,GOrder,Ctx},
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,GOrder,Ctx<:AbstractWignerContext}
    return _latex_display(io) do
        write(
            io,
            "\\begin{aligned}\nI_{\\mathrm{coll}}(k)={}&",
            _off_shell_string(collision, true),
            "\n\\end{aligned}",
        )
        return nothing
    end
end
