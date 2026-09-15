# Physics-facing displays for the staged kinetic compiler.
#
# The early diagrammatic stages are summarized by their physical representation. Once the
# collision expression exists, displays render the actual routed integrand rather than a
# description of the compiler pass that produced it.

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
    write(io, "\$\$\n")
    writer()
    write(io, "\n\$\$")
    return nothing
end

function _display_number(x::Number, latex::Bool)
    io = IOBuffer()
    if latex
        write_latex_number(io, x)
    elseif x isa Rational
        if denominator(x) == 1
            write(io, string(numerator(x)))
        else
            write(io, string(numerator(x)), "/", string(denominator(x)))
        end
    else
        show(io, x)
    end
    return String(take!(io))
end

function _display_symbol(symbol::Symbol, latex::Bool)
    latex || return string(symbol)
    io = IOBuffer()
    write_latex_symbol(io, symbol)
    return String(take!(io))
end

function _parameter_string(parameter::ParameterMonomial, latex::Bool)
    isempty(parameter.powers) && return ""
    factors = String[]
    for power in parameter.powers
        base = _display_symbol(power.name, latex)
        if power.exponent == 1
            push!(factors, base)
        elseif latex
            push!(factors, base * "^{" * string(power.exponent) * "}")
        else
            push!(factors, base * "^" * string(power.exponent))
        end
    end
    return join(factors, latex ? "\\," : " ")
end

function _momentum_label(
    basis::MomentumBasis, external::MomentumVariable, slot::Int, latex::Bool
)
    external_slot = _external_basis_index(basis, external)
    slot == external_slot && return "k"
    loop_slots = _loop_basis_indices(basis, external)
    ordinal = findfirst(==(slot), loop_slots)
    ordinal === nothing && return "q?"
    length(loop_slots) == 1 && return "q"
    return latex ? "q_{" * string(ordinal) * "}" : "q" * string(ordinal)
end

function _linear_momentum_string(
    momentum::LinearMomentum, basis::MomentumBasis, external::MomentumVariable, latex::Bool
)
    terms = String[]
    for slot in eachindex(momentum.coefficients)
        coefficient = momentum[slot]
        iszero(coefficient) && continue
        label = _momentum_label(basis, external, slot, latex)
        magnitude = abs(coefficient)
        if magnitude == 1
            body = label
        else
            multiplication = latex ? "\\," : "*"
            body = _display_number(magnitude, latex) * multiplication * label
        end
        push!(terms, coefficient < 0 ? "-" * body : body)
    end
    isempty(terms) && return "0"
    out = first(terms)
    for term in Iterators.drop(terms, 1)
        if startswith(term, "-")
            out *= " - " * term[2:end]
        else
            out *= " + " * term
        end
    end
    return out
end

function _component_string(
    component::MomentumComponent,
    basis::MomentumBasis,
    external::MomentumVariable,
    latex::Bool,
)
    momentum_text = _linear_momentum_string(component.momentum, basis, external, latex)
    nonzero_slots = findall(value -> !iszero(value), component.momentum.coefficients)
    simple =
        length(nonzero_slots) == 1 && component.momentum[first(nonzero_slots)] == 1 // 1
    axis = _display_symbol(component.axis, latex)

    if simple && latex
        slot = first(nonzero_slots)
        external_slot = _external_basis_index(basis, external)
        slot == external_slot && return "k_{" * axis * "}"

        loop_slots = _loop_basis_indices(basis, external)
        ordinal = findfirst(==(slot), loop_slots)
        ordinal === nothing && return "\\left(" * momentum_text * "\\right)_{" * axis * "}"
        if length(loop_slots) == 1
            return "q_{" * axis * "}"
        end
        return "q_{" * string(ordinal) * "," * axis * "}"
    elseif simple
        return momentum_text * "_" * axis
    elseif latex
        return "\\left(" * momentum_text * "\\right)_{" * axis * "}"
    end
    return "(" * momentum_text * ")_" * axis
end

function _grouped_factor_strings(render, factors)
    out = String[]
    i = 1
    while i <= length(factors)
        factor = factors[i]
        multiplicity = 1
        while i + multiplicity <= length(factors)
            factors[i + multiplicity] == factor || break
            multiplicity += 1
        end
        text = render(factor)
        if multiplicity > 1
            text *= "^" * string(multiplicity)
        end
        push!(out, text)
        i += multiplicity
    end
    return out
end

function _momentum_monomial_string(
    monomial::MomentumMonomial,
    basis::MomentumBasis,
    external::MomentumVariable,
    latex::Bool,
)
    factors = _grouped_factor_strings(monomial.factors) do component
        return _component_string(component, basis, external, latex)
    end
    return join(factors, latex ? "\\," : " ")
end

_family_string(family::FieldFamily, latex::Bool) = _display_symbol(name(family), latex)

function _regulator_label(shift::Int)
    iszero(shift) && return ""
    shift == 1 && return "+"
    shift == -1 && return "-"
    return shift > 0 ? "+" * string(shift) : string(shift)
end

function _regulated_family_string(family::FieldFamily, shift::Int, latex::Bool)
    family_text = _family_string(family, latex)
    regulator = _regulator_label(shift)
    isempty(regulator) && return family_text
    return family_text * "," * regulator
end

function _distribution_atom_string(
    prefix::String,
    family::FieldFamily,
    momentum::LinearMomentum,
    basis::MomentumBasis,
    external::MomentumVariable,
    latex::Bool,
)
    family_text = _family_string(family, latex)
    momentum_text = _linear_momentum_string(momentum, basis, external, latex)
    if latex
        return prefix * "_{" * family_text * "}\\!\\left(" * momentum_text * "\\right)"
    end
    return prefix * "_" * family_text * "(" * momentum_text * ")"
end

function _statistical_monomial_string(
    monomial::StatisticalMonomial,
    basis::MomentumBasis,
    external::MomentumVariable,
    latex::Bool,
)
    factors = _grouped_factor_strings(monomial.factors) do atom
        return _distribution_atom_string(
            "F", atom.family, atom.momentum, basis, external, latex
        )
    end
    return join(factors, latex ? "\\," : " ")
end

function _occupation_monomial_string(
    monomial::OccupationMonomial,
    basis::MomentumBasis,
    external::MomentumVariable,
    latex::Bool,
)
    factors = _grouped_factor_strings(monomial.factors) do atom
        return _distribution_atom_string(
            "n", atom.family, atom.momentum, basis, external, latex
        )
    end
    return join(factors, latex ? "\\," : " ")
end

function _energy_form_string(
    form::EnergyForm, basis::MomentumBasis, external::MomentumVariable, latex::Bool
)
    terms = String[]
    for (atom, coefficient) in form
        family_text = _family_string(atom.family, latex)
        momentum_text = _linear_momentum_string(atom.momentum, basis, external, latex)
        if latex
            energy =
                "\\varepsilon_{" * family_text * "}\\!\\left(" * momentum_text * "\\right)"
        else
            energy = "ε_" * family_text * "(" * momentum_text * ")"
        end
        magnitude = abs(coefficient)
        if magnitude == 1
            body = energy
        else
            multiplication = latex ? "\\," : "*"
            body = _display_number(magnitude, latex) * multiplication * energy
        end
        push!(terms, coefficient < 0 ? "-" * body : body)
    end
    isempty(terms) && return "0"
    out = first(terms)
    for term in Iterators.drop(terms, 1)
        if startswith(term, "-")
            out *= " - " * term[2:end]
        else
            out *= " + " * term
        end
    end
    return out
end

function _frequency_support_string(
    support::FrequencySupport, basis::MomentumBasis, external::MomentumVariable, latex::Bool
)
    factors = String[]
    for shell in support.shells
        energy = _energy_form_string(shell.energy, basis, external, latex)
        if latex
            push!(factors, "\\delta\\!\\left(" * energy * "\\right)")
        else
            push!(factors, "δ(" * energy * ")")
        end
    end
    for pv in support.principal_values
        energy = _energy_form_string(pv.energy, basis, external, latex)
        if latex
            push!(factors, "\\operatorname{PV}\\!\\left(\\frac{1}{" * energy * "}\\right)")
        else
            push!(factors, "PV(1/(" * energy * "))")
        end
    end
    return join(factors, latex ? "\\," : " ")
end

function _measure_string(basis::MomentumBasis, external::MomentumVariable, latex::Bool)
    loops = _loop_basis_indices(basis, external)
    labels = [_momentum_label(basis, external, slot, latex) for slot in loops]
    if latex
        return join(["\\int_{" * label * "}" for label in labels])
    end
    return join(["∫_" * label for label in labels], " ")
end

function _term_string(coefficient::Number, factors::Vector{String}, latex::Bool)
    separator = latex ? "\\," : " "
    nonempty_factors = [factor for factor in factors if !isempty(factor)]
    body = join(nonempty_factors, separator)
    isempty(body) && return _display_number(coefficient, latex)
    isone(coefficient) && return body
    isone(-coefficient) && return "-" * body
    return _display_number(coefficient, latex) * separator * body
end

function _sum_string(terms::Vector{String}, latex::Bool)
    isempty(terms) && return "0"
    separator = latex ? "\\\\\n&{}" : "\n    "
    out = first(terms)
    for term in Iterators.drop(terms, 1)
        if startswith(term, "-")
            out *= separator * "- " * term[2:end]
        else
            out *= separator * "+ " * term
        end
    end
    return out
end

function _polynomial_monomial_string(
    monomial::StatisticalMonomial,
    basis::MomentumBasis,
    external::MomentumVariable,
    latex::Bool,
)
    return _statistical_monomial_string(monomial, basis, external, latex)
end

function _polynomial_monomial_string(
    monomial::OccupationMonomial,
    basis::MomentumBasis,
    external::MomentumVariable,
    latex::Bool,
)
    return _occupation_monomial_string(monomial, basis, external, latex)
end

function _physical_collision_string(terms, latex::Bool)
    rendered = String[]
    for (sector, polynomial) in terms
        basis = momentum_basis(sector)
        external = external_wigner_momentum(sector)
        parameter_text = _parameter_string(parameters(sector), latex)
        measure_text = _measure_string(basis, external, latex)
        support_text = _frequency_support_string(
            frequency_support(sector), basis, external, latex
        )
        for (kinematic, kinematic_coefficient) in kinematic_factor(sector)
            kinematic_text = _momentum_monomial_string(kinematic, basis, external, latex)
            for (monomial, polynomial_coefficient) in polynomial
                coefficient = kinematic_coefficient * polynomial_coefficient
                iszero(coefficient) && continue
                polynomial_text = _polynomial_monomial_string(
                    monomial, basis, external, latex
                )
                factors = String[
                    parameter_text,
                    measure_text,
                    kinematic_text,
                    polynomial_text,
                    support_text,
                ]
                push!(rendered, _term_string(coefficient, factors, latex))
            end
        end
    end
    return _sum_string(rendered, latex)
end

function _sd_line_string(
    line::KineticLine,
    kind::SpectralDispersiveKind,
    basis::MomentumBasis,
    external::MomentumVariable,
    latex::Bool,
)
    prefix = kind === CollisionSpectral ? "A" : "D"
    shift = Int(regularisation_shift(line))
    family = _regulated_family_string(line.family, shift, latex)
    momentum_text = _linear_momentum_string(momentum(line), basis, external, latex)
    if latex
        return prefix * "_{" * family * "}\\!\\left(" * momentum_text * "\\right)"
    end
    return prefix * "_" * family * "(" * momentum_text * ")"
end

function _sd_expression_terms(
    collision::SpectralDispersiveCollision,
    expression::SpectralDispersiveExpression,
    include_external::Bool,
    latex::Bool,
)
    rendered = String[]
    parameter_text = _parameter_string(parameters(collision), latex)
    for (term, source_coefficient) in expression
        basis = momentum_basis(term)
        external = external_wigner_momentum(term)
        measure_text = _measure_string(basis, external, latex)
        lines = kinetic_lines(term.carrier)
        kinds = spectral_dispersive_kinds(term)
        line_factors = String[]
        if include_external
            external_slot = _external_basis_index(basis, external)
            external_momentum = basis_momentum(basis, external_slot)
            push!(
                line_factors,
                _distribution_atom_string(
                    "F", target_family(collision), external_momentum, basis, external, latex
                ),
            )
        end
        for i in eachindex(lines)
            line = lines[i]
            if statistical_weight(line) === DistributionWeight
                push!(
                    line_factors,
                    _distribution_atom_string(
                        "F", line.family, momentum(line), basis, external, latex
                    ),
                )
            end
            push!(line_factors, _sd_line_string(line, kinds[i], basis, external, latex))
        end
        for (kinematic, kinematic_coefficient) in kinematic_factor(term)
            coefficient = source_coefficient * kinematic_coefficient
            iszero(coefficient) && continue
            kinematic_text = _momentum_monomial_string(kinematic, basis, external, latex)
            factors = vcat(
                String[parameter_text, measure_text, kinematic_text], line_factors
            )
            push!(rendered, _term_string(coefficient, factors, latex))
        end
    end
    return rendered
end

function _spectral_dispersive_string(collision::SpectralDispersiveCollision, latex::Bool)
    rendered = String[]
    append!(
        rendered, _sd_expression_terms(collision, collision_offset(collision), false, latex)
    )
    append!(
        rendered,
        _sd_expression_terms(
            collision, collision_distribution_coefficient(collision), true, latex
        ),
    )
    return _sum_string(rendered, latex)
end

function _blocked_summary(result::ReducedFrequencyCollision, latex::Bool)
    blocked = reduced_blocked_terms(result)
    unresolved = reduced_trotter_terms(result)
    if latex
        return "N_{\\mathrm{blocked}}=" *
               string(length(blocked)) *
               ",\\qquad N_{\\mathrm{Trotter}}=" *
               string(length(unresolved))
    end
    return "causal blockers: " *
           string(length(blocked)) *
           "\n  unresolved Trotter states: " *
           string(length(unresolved))
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
    write(io, "\n  I_coll(k) = ", _spectral_dispersive_string(collision, false))
    return nothing
end

function Base.show(io::IO, ::MIME"text/latex", collision::SpectralDispersiveCollision)
    return _latex_display(io) do
        write(
            io,
            "\\begin{aligned}\nI_{\\mathrm{coll}}(k)={}&",
            _spectral_dispersive_string(collision, true),
            "\n\\end{aligned}",
        )
        return nothing
    end
end
function Base.show(io::IO, ::MIME"text/plain", result::ReducedFrequencyCollision)
    _show_stage_plain(io, "Reduced causal-frequency collision", result)
    write(
        io,
        "\n  I_reg(k) = ",
        _physical_collision_string(reduced_regular_terms(result), false),
    )
    write(io, "\n  ", _blocked_summary(result, false))
    return nothing
end

function Base.show(io::IO, ::MIME"text/latex", result::ReducedFrequencyCollision)
    return _latex_display(io) do
        write(
            io,
            "\\begin{aligned}\nI_{\\mathrm{reg}}(k)={}&",
            _physical_collision_string(reduced_regular_terms(result), true),
            "\\\\\n&",
            _blocked_summary(result, true),
            "\n\\end{aligned}",
        )
        return nothing
    end
end

function Base.show(io::IO, ::MIME"text/plain", result::OccupationReducedExpression)
    _show_stage_plain(io, "Occupation-reduced regular collision", result)
    write(
        io,
        "\n  C_n^reg(k) = ",
        _physical_collision_string(occupation_reduced_terms(result), false),
    )
    return nothing
end

function Base.show(io::IO, ::MIME"text/latex", result::OccupationReducedExpression)
    return _latex_display(io) do
        write(
            io,
            "\\begin{aligned}\nC_n^{\\mathrm{reg}}(k)={}&",
            _physical_collision_string(occupation_reduced_terms(result), true),
            "\n\\end{aligned}",
        )
        return nothing
    end
end

function Base.show(io::IO, ::MIME"text/plain", result::LoopQuotientedExpression)
    _show_stage_plain(io, "Loop-quotiented regular collision", result)
    write(
        io,
        "\n  C_n^quot(k) = ",
        _physical_collision_string(loop_quotient_terms(result), false),
    )
    return nothing
end

function Base.show(io::IO, ::MIME"text/latex", result::LoopQuotientedExpression)
    return _latex_display(io) do
        write(
            io,
            "\\begin{aligned}\nC_n^{\\mathrm{quot}}(k)={}&",
            _physical_collision_string(loop_quotient_terms(result), true),
            "\n\\end{aligned}",
        )
        return nothing
    end
end

function Base.show(io::IO, ::MIME"text/plain", kernel::CollisionKernel)
    _show_stage_plain(io, "Collision kernel", kernel)
    write(
        io,
        "\n  C_n(k) = ",
        _physical_collision_string(collision_kernel_terms(kernel), false),
    )
    return nothing
end

function Base.show(io::IO, ::MIME"text/latex", kernel::CollisionKernel)
    return _latex_display(io) do
        write(
            io,
            "\\begin{aligned}\nC_n(k)={}&",
            _physical_collision_string(collision_kernel_terms(kernel), true),
            "\n\\end{aligned}",
        )
        return nothing
    end
end