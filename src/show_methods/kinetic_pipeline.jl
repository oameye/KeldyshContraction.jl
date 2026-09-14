# Physics-facing displays for the staged kinetic compiler.
#
# Stage displays must show the represented physical expression, not merely restate the
# transformation that produced it. The small metadata summaries below are retained for the
# earlier diagrammatic stages; the collision stages render their actual routed integrands.

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

function _write_plain_parameter(io::IO, parameter::ParameterMonomial)
    if isempty(parameter.powers)
        write(io, "1")
        return nothing
    end
    for (i, power) in enumerate(parameter.powers)
        i > 1 && write(io, " ")
        write(io, string(power.name))
        power.exponent == 1 || write(io, "^", string(power.exponent))
    end
    return nothing
end

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

function _write_plain_number(io::IO, x::Rational)
    if denominator(x) == 1
        write(io, string(numerator(x)))
    else
        write(io, string(numerator(x)), "/", string(denominator(x)))
    end
    return nothing
end

function _write_plain_number(io::IO, x::Real)
    show(io, x)
    return nothing
end

function _write_plain_number(io::IO, x::Complex)
    re = real(x)
    im_part = imag(x)
    if iszero(im_part)
        _write_plain_number(io, re)
    elseif iszero(re)
        if isone(im_part)
            write(io, "i")
        elseif isone(-im_part)
            write(io, "-i")
        else
            _write_plain_number(io, im_part)
            write(io, "i")
        end
    else
        write(io, "(")
        _write_plain_number(io, re)
        if im_part < 0
            write(io, " - ")
            _write_plain_number(io, -im_part)
        else
            write(io, " + ")
            _write_plain_number(io, im_part)
        end
        write(io, "i)")
    end
    return nothing
end

function _write_plain_number(io::IO, x::Number)
    show(io, x)
    return nothing
end

_plain_negative_real(x::Real) = x < 0
function _plain_negative_real(x::Complex)
    return iszero(imag(x)) && real(x) < 0
end
_plain_negative_real(::Number) = false

_parameter_isone(parameter::ParameterMonomial) = isempty(parameter.powers)

function _basis_variable_slot(
    basis::MomentumBasis, external::MomentumVariable, variable::MomentumVariable
)
    for i in eachindex(basis.variables)
        basis[i] == variable && return i
    end
    throw(ArgumentError("momentum variable is not present in the displayed basis"))
end

function _loop_ordinal(basis::MomentumBasis, external::MomentumVariable, slot::Int)
    external_slot = _external_basis_index(basis, external)
    slot == external_slot && return 0
    ordinal = 0
    for i in eachindex(basis.variables)
        i == external_slot && continue
        ordinal += 1
        i == slot && return ordinal
    end
    throw(ArgumentError("momentum slot is not present in the displayed loop basis"))
end

function _write_plain_basis_variable(
    io::IO, basis::MomentumBasis, external::MomentumVariable, slot::Int
)
    external_slot = _external_basis_index(basis, external)
    if slot == external_slot
        write(io, "k")
        return nothing
    end
    nloops = length(basis) - 1
    ordinal = _loop_ordinal(basis, external, slot)
    write(io, nloops == 1 ? "q" : "q" * string(ordinal))
    return nothing
end

function _write_latex_basis_variable(
    io::IO, basis::MomentumBasis, external::MomentumVariable, slot::Int
)
    external_slot = _external_basis_index(basis, external)
    if slot == external_slot
        write(io, "k")
        return nothing
    end
    nloops = length(basis) - 1
    ordinal = _loop_ordinal(basis, external, slot)
    if nloops == 1
        write(io, "q")
    else
        write(io, "q_{", string(ordinal), "}")
    end
    return nothing
end

function _write_plain_linear_momentum(
    io::IO, momentum::LinearMomentum, basis::MomentumBasis, external::MomentumVariable
)
    first_term = true
    for slot in eachindex(momentum.coefficients)
        coefficient = momentum[slot]
        iszero(coefficient) && continue
        negative = coefficient < 0
        magnitude = abs(coefficient)
        if first_term
            negative && write(io, "-")
        else
            write(io, negative ? " - " : " + ")
        end
        if magnitude != 1
            _write_plain_number(io, magnitude)
            write(io, "*")
        end
        _write_plain_basis_variable(io, basis, external, slot)
        first_term = false
    end
    first_term && write(io, "0")
    return nothing
end

function _write_latex_linear_momentum(
    io::IO, momentum::LinearMomentum, basis::MomentumBasis, external::MomentumVariable
)
    first_term = true
    for slot in eachindex(momentum.coefficients)
        coefficient = momentum[slot]
        iszero(coefficient) && continue
        negative = coefficient < 0
        magnitude = abs(coefficient)
        if first_term
            negative && write(io, "-")
        else
            write(io, negative ? "-" : "+")
        end
        if magnitude != 1
            write_latex_number(io, magnitude)
            write(io, "\\,")
        end
        _write_latex_basis_variable(io, basis, external, slot)
        first_term = false
    end
    first_term && write(io, "0")
    return nothing
end

function _simple_unit_momentum_slot(momentum::LinearMomentum)
    found = 0
    for slot in eachindex(momentum.coefficients)
        coefficient = momentum[slot]
        iszero(coefficient) && continue
        coefficient == 1 || return 0
        iszero(found) || return 0
        found = slot
    end
    return found
end

function _write_plain_momentum_component(
    io::IO, component::MomentumComponent, basis::MomentumBasis, external::MomentumVariable
)
    slot = _simple_unit_momentum_slot(component.momentum)
    if iszero(slot)
        write(io, "(")
        _write_plain_linear_momentum(io, component.momentum, basis, external)
        write(io, ")")
    else
        _write_plain_basis_variable(io, basis, external, slot)
    end
    write(io, "_", string(component.axis))
    return nothing
end

function _write_latex_momentum_component(
    io::IO, component::MomentumComponent, basis::MomentumBasis, external::MomentumVariable
)
    slot = _simple_unit_momentum_slot(component.momentum)
    if iszero(slot)
        write(io, "\\left(")
        _write_latex_linear_momentum(io, component.momentum, basis, external)
        write(io, "\\right)")
    else
        _write_latex_basis_variable(io, basis, external, slot)
    end
    write(io, "_{")
    write_latex_symbol(io, component.axis)
    write(io, "}")
    return nothing
end

function _write_plain_momentum_monomial(
    io::IO, monomial::MomentumMonomial, basis::MomentumBasis, external::MomentumVariable
)
    i = 1
    while i <= length(monomial.factors)
        component = monomial.factors[i]
        multiplicity = 1
        while i + multiplicity <= length(monomial.factors) &&
            monomial.factors[i + multiplicity] == component
            multiplicity += 1
        end
        i > 1 && write(io, " ")
        _write_plain_momentum_component(io, component, basis, external)
        multiplicity > 1 && write(io, "^", string(multiplicity))
        i += multiplicity
    end
    return nothing
end

function _write_latex_momentum_monomial(
    io::IO, monomial::MomentumMonomial, basis::MomentumBasis, external::MomentumVariable
)
    i = 1
    while i <= length(monomial.factors)
        component = monomial.factors[i]
        multiplicity = 1
        while i + multiplicity <= length(monomial.factors) &&
            monomial.factors[i + multiplicity] == component
            multiplicity += 1
        end
        i > 1 && write(io, "\\,")
        _write_latex_momentum_component(io, component, basis, external)
        multiplicity > 1 && write(io, "^{", string(multiplicity), "}")
        i += multiplicity
    end
    return nothing
end

function _write_plain_family(io::IO, family::FieldFamily)
    write(io, string(name(family)))
    return nothing
end

function _write_latex_family(io::IO, family::FieldFamily)
    write_latex_symbol(io, name(family))
    return nothing
end

function _write_plain_statistical_atom(
    io::IO, atom::StatisticalAtom, basis::MomentumBasis, external::MomentumVariable
)
    write(io, "F_")
    _write_plain_family(io, atom.family)
    write(io, "(")
    _write_plain_linear_momentum(io, atom.momentum, basis, external)
    write(io, ")")
    return nothing
end

function _write_latex_statistical_atom(
    io::IO, atom::StatisticalAtom, basis::MomentumBasis, external::MomentumVariable
)
    write(io, "F_{")
    _write_latex_family(io, atom.family)
    write(io, "}\\!\\left(")
    _write_latex_linear_momentum(io, atom.momentum, basis, external)
    write(io, "\\right)")
    return nothing
end

function _write_plain_occupation_atom(
    io::IO, atom::OccupationAtom, basis::MomentumBasis, external::MomentumVariable
)
    write(io, "n_")
    _write_plain_family(io, atom.family)
    write(io, "(")
    _write_plain_linear_momentum(io, atom.momentum, basis, external)
    write(io, ")")
    return nothing
end

function _write_latex_occupation_atom(
    io::IO, atom::OccupationAtom, basis::MomentumBasis, external::MomentumVariable
)
    write(io, "n_{")
    _write_latex_family(io, atom.family)
    write(io, "}\\!\\left(")
    _write_latex_linear_momentum(io, atom.momentum, basis, external)
    write(io, "\\right)")
    return nothing
end

function _write_plain_statistical_monomial(
    io::IO, monomial::StatisticalMonomial, basis::MomentumBasis, external::MomentumVariable
)
    i = 1
    while i <= length(monomial.factors)
        atom = monomial.factors[i]
        multiplicity = 1
        while i + multiplicity <= length(monomial.factors) &&
            monomial.factors[i + multiplicity] == atom
            multiplicity += 1
        end
        i > 1 && write(io, " ")
        _write_plain_statistical_atom(io, atom, basis, external)
        multiplicity > 1 && write(io, "^", string(multiplicity))
        i += multiplicity
    end
    return nothing
end

function _write_latex_statistical_monomial(
    io::IO, monomial::StatisticalMonomial, basis::MomentumBasis, external::MomentumVariable
)
    i = 1
    while i <= length(monomial.factors)
        atom = monomial.factors[i]
        multiplicity = 1
        while i + multiplicity <= length(monomial.factors) &&
            monomial.factors[i + multiplicity] == atom
            multiplicity += 1
        end
        i > 1 && write(io, "\\,")
        _write_latex_statistical_atom(io, atom, basis, external)
        multiplicity > 1 && write(io, "^{", string(multiplicity), "}")
        i += multiplicity
    end
    return nothing
end

function _write_plain_occupation_monomial(
    io::IO, monomial::OccupationMonomial, basis::MomentumBasis, external::MomentumVariable
)
    i = 1
    while i <= length(monomial.factors)
        atom = monomial.factors[i]
        multiplicity = 1
        while i + multiplicity <= length(monomial.factors) &&
            monomial.factors[i + multiplicity] == atom
            multiplicity += 1
        end
        i > 1 && write(io, " ")
        _write_plain_occupation_atom(io, atom, basis, external)
        multiplicity > 1 && write(io, "^", string(multiplicity))
        i += multiplicity
    end
    return nothing
end

function _write_latex_occupation_monomial(
    io::IO, monomial::OccupationMonomial, basis::MomentumBasis, external::MomentumVariable
)
    i = 1
    while i <= length(monomial.factors)
        atom = monomial.factors[i]
        multiplicity = 1
        while i + multiplicity <= length(monomial.factors) &&
            monomial.factors[i + multiplicity] == atom
            multiplicity += 1
        end
        i > 1 && write(io, "\\,")
        _write_latex_occupation_atom(io, atom, basis, external)
        multiplicity > 1 && write(io, "^{", string(multiplicity), "}")
        i += multiplicity
    end
    return nothing
end

function _write_plain_dispersion_atom(
    io::IO, atom::DispersionAtom, basis::MomentumBasis, external::MomentumVariable
)
    write(io, "ε_")
    _write_plain_family(io, atom.family)
    write(io, "(")
    _write_plain_linear_momentum(io, atom.momentum, basis, external)
    write(io, ")")
    return nothing
end

function _write_latex_dispersion_atom(
    io::IO, atom::DispersionAtom, basis::MomentumBasis, external::MomentumVariable
)
    write(io, "\\varepsilon_{")
    _write_latex_family(io, atom.family)
    write(io, "}\\!\\left(")
    _write_latex_linear_momentum(io, atom.momentum, basis, external)
    write(io, "\\right)")
    return nothing
end

function _write_plain_energy_form(
    io::IO, form::EnergyForm, basis::MomentumBasis, external::MomentumVariable
)
    first_term = true
    for (atom, coefficient) in form
        negative = coefficient < 0
        magnitude = abs(coefficient)
        if first_term
            negative && write(io, "-")
        else
            write(io, negative ? " - " : " + ")
        end
        if magnitude != 1
            _write_plain_number(io, magnitude)
            write(io, "*")
        end
        _write_plain_dispersion_atom(io, atom, basis, external)
        first_term = false
    end
    first_term && write(io, "0")
    return nothing
end

function _write_latex_energy_form(
    io::IO, form::EnergyForm, basis::MomentumBasis, external::MomentumVariable
)
    first_term = true
    for (atom, coefficient) in form
        negative = coefficient < 0
        magnitude = abs(coefficient)
        if first_term
            negative && write(io, "-")
        else
            write(io, negative ? "-" : "+")
        end
        if magnitude != 1
            write_latex_number(io, magnitude)
            write(io, "\\,")
        end
        _write_latex_dispersion_atom(io, atom, basis, external)
        first_term = false
    end
    first_term && write(io, "0")
    return nothing
end

_support_isempty(support::FrequencySupport) =
    isempty(support.shells) && isempty(support.principal_values)

function _write_plain_frequency_support(
    io::IO, support::FrequencySupport, basis::MomentumBasis, external::MomentumVariable
)
    wrote = false
    for shell in support.shells
        wrote && write(io, " ")
        write(io, "δ(")
        _write_plain_energy_form(io, shell.energy, basis, external)
        write(io, ")")
        wrote = true
    end
    for principal_value in support.principal_values
        wrote && write(io, " ")
        write(io, "PV(1/(")
        _write_plain_energy_form(io, principal_value.energy, basis, external)
        write(io, "))")
        wrote = true
    end
    return nothing
end

function _write_latex_frequency_support(
    io::IO, support::FrequencySupport, basis::MomentumBasis, external::MomentumVariable
)
    wrote = false
    for shell in support.shells
        wrote && write(io, "\\,")
        write(io, "\\delta\\!\\left(")
        _write_latex_energy_form(io, shell.energy, basis, external)
        write(io, "\\right)")
        wrote = true
    end
    for principal_value in support.principal_values
        wrote && write(io, "\\,")
        write(io, "\\operatorname{PV}\\!\\left(\\frac{1}{")
        _write_latex_energy_form(io, principal_value.energy, basis, external)
        write(io, "}\\right)")
        wrote = true
    end
    return nothing
end

function _write_plain_measure(
    io::IO, basis::MomentumBasis, external::MomentumVariable
)
    loops = _loop_basis_indices(basis, external)
    for (i, slot) in enumerate(loops)
        i > 1 && write(io, " ")
        write(io, "∫_")
        _write_plain_basis_variable(io, basis, external, slot)
    end
    return !isempty(loops)
end

function _write_latex_measure(
    io::IO, basis::MomentumBasis, external::MomentumVariable
)
    loops = _loop_basis_indices(basis, external)
    for slot in loops
        write(io, "\\int_{")
        _write_latex_basis_variable(io, basis, external, slot)
        write(io, "}")
    end
    return !isempty(loops)
end

function _write_plain_contribution(
    io::IO,
    coefficient::Number,
    parameter::ParameterMonomial,
    basis::MomentumBasis,
    external::MomentumVariable,
    kinematic::MomentumMonomial,
    statistical::Union{StatisticalMonomial,OccupationMonomial},
    support::Union{FrequencySupport,Nothing},
)
    wrote = false
    has_non_numeric =
        !_parameter_isone(parameter) ||
        length(basis) > 1 ||
        !isempty(kinematic) ||
        !isempty(statistical) ||
        (support !== nothing && !_support_isempty(support))
    if !isone(coefficient) || !has_non_numeric
        _write_plain_number(io, coefficient)
        wrote = true
    end
    if !_parameter_isone(parameter)
        wrote && write(io, " ")
        _write_plain_parameter(io, parameter)
        wrote = true
    end
    if length(basis) > 1
        wrote && write(io, " ")
        _write_plain_measure(io, basis, external)
        wrote = true
    end
    if !isempty(kinematic)
        wrote && write(io, " ")
        _write_plain_momentum_monomial(io, kinematic, basis, external)
        wrote = true
    end
    if !isempty(statistical)
        wrote && write(io, " ")
        if statistical isa StatisticalMonomial
            _write_plain_statistical_monomial(io, statistical, basis, external)
        else
            _write_plain_occupation_monomial(io, statistical, basis, external)
        end
        wrote = true
    end
    if support !== nothing && !_support_isempty(support)
        wrote && write(io, " ")
        _write_plain_frequency_support(io, support, basis, external)
    end
    return nothing
end

function _write_latex_contribution(
    io::IO,
    coefficient::Number,
    parameter::ParameterMonomial,
    basis::MomentumBasis,
    external::MomentumVariable,
    kinematic::MomentumMonomial,
    statistical::Union{StatisticalMonomial,OccupationMonomial},
    support::Union{FrequencySupport,Nothing},
)
    wrote = false
    has_non_numeric =
        !_parameter_isone(parameter) ||
        length(basis) > 1 ||
        !isempty(kinematic) ||
        !isempty(statistical) ||
        (support !== nothing && !_support_isempty(support))
    if !isone(coefficient) || !has_non_numeric
        write_latex_number(io, coefficient)
        wrote = true
    end
    if !_parameter_isone(parameter)
        wrote && write(io, "\\,")
        _write_latex_parameter(io, parameter)
        wrote = true
    end
    if length(basis) > 1
        wrote && write(io, "\\,")
        _write_latex_measure(io, basis, external)
        wrote = true
    end
    if !isempty(kinematic)
        wrote && write(io, "\\,")
        _write_latex_momentum_monomial(io, kinematic, basis, external)
        wrote = true
    end
    if !isempty(statistical)
        wrote && write(io, "\\,")
        if statistical isa StatisticalMonomial
            _write_latex_statistical_monomial(io, statistical, basis, external)
        else
            _write_latex_occupation_monomial(io, statistical, basis, external)
        end
        wrote = true
    end
    if support !== nothing && !_support_isempty(support)
        wrote && write(io, "\\,")
        _write_latex_frequency_support(io, support, basis, external)
    end
    return nothing
end

function _write_plain_signed_prefix(io::IO, coefficient::Number, first_term::Bool)
    negative = _plain_negative_real(coefficient)
    magnitude = negative ? -coefficient : coefficient
    if first_term
        negative && write(io, "-")
    else
        write(io, negative ? "\n    - " : "\n    + ")
    end
    return magnitude
end

function _write_latex_signed_prefix(io::IO, coefficient::Number, first_term::Bool)
    negative = _latex_negative_real(coefficient)
    magnitude = negative ? -coefficient : coefficient
    if first_term
        negative && write(io, "-")
    else
        write(io, negative ? "\\\\\n&{}-" : "\\\\\n&{}+")
    end
    return magnitude
end

function _write_plain_statistical_collision(io::IO, terms)
    first_term = true
    for (sector, polynomial) in terms
        basis = momentum_basis(sector)
        external = external_wigner_momentum(sector)
        for (kinematic, kinematic_coefficient) in kinematic_factor(sector)
            for (statistical, statistical_coefficient) in polynomial
                coefficient = kinematic_coefficient * statistical_coefficient
                iszero(coefficient) && continue
                magnitude = _write_plain_signed_prefix(io, coefficient, first_term)
                _write_plain_contribution(
                    io,
                    magnitude,
                    parameters(sector),
                    basis,
                    external,
                    kinematic,
                    statistical,
                    frequency_support(sector),
                )
                first_term = false
            end
        end
    end
    first_term && write(io, "0")
    return nothing
end

function _write_latex_statistical_collision(io::IO, terms)
    first_term = true
    for (sector, polynomial) in terms
        basis = momentum_basis(sector)
        external = external_wigner_momentum(sector)
        for (kinematic, kinematic_coefficient) in kinematic_factor(sector)
            for (statistical, statistical_coefficient) in polynomial
                coefficient = kinematic_coefficient * statistical_coefficient
                iszero(coefficient) && continue
                magnitude = _write_latex_signed_prefix(io, coefficient, first_term)
                _write_latex_contribution(
                    io,
                    magnitude,
                    parameters(sector),
                    basis,
                    external,
                    kinematic,
                    statistical,
                    frequency_support(sector),
                )
                first_term = false
            end
        end
    end
    first_term && write(io, "0")
    return nothing
end

function _write_plain_occupation_collision(io::IO, terms)
    first_term = true
    for (sector, polynomial) in terms
        basis = momentum_basis(sector)
        external = external_wigner_momentum(sector)
        for (kinematic, kinematic_coefficient) in kinematic_factor(sector)
            for (occupation, occupation_coefficient) in polynomial
                coefficient = kinematic_coefficient * occupation_coefficient
                iszero(coefficient) && continue
                magnitude = _write_plain_signed_prefix(io, coefficient, first_term)
                _write_plain_contribution(
                    io,
                    magnitude,
                    parameters(sector),
                    basis,
                    external,
                    kinematic,
                    occupation,
                    frequency_support(sector),
                )
                first_term = false
            end
        end
    end
    first_term && write(io, "0")
    return nothing
end

function _write_latex_occupation_collision(io::IO, terms)
    first_term = true
    for (sector, polynomial) in terms
        basis = momentum_basis(sector)
        external = external_wigner_momentum(sector)
        for (kinematic, kinematic_coefficient) in kinematic_factor(sector)
            for (occupation, occupation_coefficient) in polynomial
                coefficient = kinematic_coefficient * occupation_coefficient
                iszero(coefficient) && continue
                magnitude = _write_latex_signed_prefix(io, coefficient, first_term)
                _write_latex_contribution(
                    io,
                    magnitude,
                    parameters(sector),
                    basis,
                    external,
                    kinematic,
                    occupation,
                    frequency_support(sector),
                )
                first_term = false
            end
        end
    end
    first_term && write(io, "0")
    return nothing
end

function _write_plain_sd_line(
    io::IO,
    line::KineticLine,
    kind::SpectralDispersiveKind,
    basis::MomentumBasis,
    external::MomentumVariable,
)
    write(io, kind === CollisionSpectral ? "A_" : "D_")
    _write_plain_family(io, line.family)
    shift = Int(regularisation_shift(line))
    if !iszero(shift)
        write(io, "[Δt=", shift > 0 ? "+" : "", string(shift), "·0+]")
    end
    write(io, "(")
    _write_plain_linear_momentum(io, momentum(line), basis, external)
    write(io, ")")
    return nothing
end

function _write_latex_sd_line(
    io::IO,
    line::KineticLine,
    kind::SpectralDispersiveKind,
    basis::MomentumBasis,
    external::MomentumVariable,
)
    write(io, kind === CollisionSpectral ? "A_{" : "D_{")
    _write_latex_family(io, line.family)
    write(io, "}")
    shift = Int(regularisation_shift(line))
    if !iszero(shift)
        write(io, "^{[\\Delta t=", shift > 0 ? "+" : "", string(shift), "\\,0^+]}")
    end
    write(io, "\\!\\left(")
    _write_latex_linear_momentum(io, momentum(line), basis, external)
    write(io, "\\right)")
    return nothing
end

function _write_plain_sd_contribution(
    io::IO,
    coefficient::Number,
    parameter::ParameterMonomial,
    term::SpectralDispersiveTerm,
    kinematic::MomentumMonomial,
    include_external::Bool,
    target::FieldFamily,
)
    basis = momentum_basis(term)
    external = external_wigner_momentum(term)
    has_non_numeric =
        !_parameter_isone(parameter) || length(basis) > 1 || !isempty(kinematic) ||
        include_external || !isempty(kinetic_lines(term.carrier))
    wrote = false
    if !isone(coefficient) || !has_non_numeric
        _write_plain_number(io, coefficient)
        wrote = true
    end
    if !_parameter_isone(parameter)
        wrote && write(io, " ")
        _write_plain_parameter(io, parameter)
        wrote = true
    end
    if length(basis) > 1
        wrote && write(io, " ")
        _write_plain_measure(io, basis, external)
        wrote = true
    end
    if !isempty(kinematic)
        wrote && write(io, " ")
        _write_plain_momentum_monomial(io, kinematic, basis, external)
        wrote = true
    end
    if include_external
        wrote && write(io, " ")
        external_slot = _external_basis_index(basis, external)
        atom = StatisticalAtom(target, basis_momentum(basis, external_slot))
        _write_plain_statistical_atom(io, atom, basis, external)
        wrote = true
    end
    lines = kinetic_lines(term.carrier)
    kinds = spectral_dispersive_kinds(term)
    for i in eachindex(lines)
        line = lines[i]
        if statistical_weight(line) === DistributionWeight
            wrote && write(io, " ")
            atom = StatisticalAtom(line.family, momentum(line))
            _write_plain_statistical_atom(io, atom, basis, external)
            wrote = true
        end
        wrote && write(io, " ")
        _write_plain_sd_line(io, line, kinds[i], basis, external)
        wrote = true
    end
    return nothing
end

function _write_latex_sd_contribution(
    io::IO,
    coefficient::Number,
    parameter::ParameterMonomial,
    term::SpectralDispersiveTerm,
    kinematic::MomentumMonomial,
    include_external::Bool,
    target::FieldFamily,
)
    basis = momentum_basis(term)
    external = external_wigner_momentum(term)
    has_non_numeric =
        !_parameter_isone(parameter) || length(basis) > 1 || !isempty(kinematic) ||
        include_external || !isempty(kinetic_lines(term.carrier))
    wrote = false
    if !isone(coefficient) || !has_non_numeric
        write_latex_number(io, coefficient)
        wrote = true
    end
    if !_parameter_isone(parameter)
        wrote && write(io, "\\,")
        _write_latex_parameter(io, parameter)
        wrote = true
    end
    if length(basis) > 1
        wrote && write(io, "\\,")
        _write_latex_measure(io, basis, external)
        wrote = true
    end
    if !isempty(kinematic)
        wrote && write(io, "\\,")
        _write_latex_momentum_monomial(io, kinematic, basis, external)
        wrote = true
    end
    if include_external
        wrote && write(io, "\\,")
        external_slot = _external_basis_index(basis, external)
        atom = StatisticalAtom(target, basis_momentum(basis, external_slot))
        _write_latex_statistical_atom(io, atom, basis, external)
        wrote = true
    end
    lines = kinetic_lines(term.carrier)
    kinds = spectral_dispersive_kinds(term)
    for i in eachindex(lines)
        line = lines[i]
        if statistical_weight(line) === DistributionWeight
            wrote && write(io, "\\,")
            atom = StatisticalAtom(line.family, momentum(line))
            _write_latex_statistical_atom(io, atom, basis, external)
            wrote = true
        end
        wrote && write(io, "\\,")
        _write_latex_sd_line(io, line, kinds[i], basis, external)
        wrote = true
    end
    return nothing
end

function _write_plain_sd_expression(io::IO, collision::SpectralDispersiveCollision)
    first_term = true
    parameter = parameters(collision)
    target = target_family(collision)
    for (expression, include_external) in (
        (collision_offset(collision), false),
        (collision_distribution_coefficient(collision), true),
    )
        for (term, source_coefficient) in expression
            for (kinematic, kinematic_coefficient) in kinematic_factor(term)
                coefficient = source_coefficient * kinematic_coefficient
                iszero(coefficient) && continue
                magnitude = _write_plain_signed_prefix(io, coefficient, first_term)
                _write_plain_sd_contribution(
                    io, magnitude, parameter, term, kinematic, include_external, target
                )
                first_term = false
            end
        end
    end
    first_term && write(io, "0")
    return nothing
end

function _write_latex_sd_expression(io::IO, collision::SpectralDispersiveCollision)
    first_term = true
    parameter = parameters(collision)
    target = target_family(collision)
    for (expression, include_external) in (
        (collision_offset(collision), false),
        (collision_distribution_coefficient(collision), true),
    )
        for (term, source_coefficient) in expression
            for (kinematic, kinematic_coefficient) in kinematic_factor(term)
                coefficient = source_coefficient * kinematic_coefficient
                iszero(coefficient) && continue
                magnitude = _write_latex_signed_prefix(io, coefficient, first_term)
                _write_latex_sd_contribution(
                    io, magnitude, parameter, term, kinematic, include_external, target
                )
                first_term = false
            end
        end
    end
    first_term && write(io, "0")
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
    write(io, "\n  I_coll(k) = ")
    _write_plain_sd_expression(io, collision)
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", collision::SpectralDispersiveCollision)
    return _latex_display(io) do
        write(io, "\\begin{aligned}\nI_{\\mathrm{coll}}(k)={}&")
        _write_latex_sd_expression(io, collision)
        write(io, "\n\\end{aligned}")
        return nothing
    end
end

function Base.show(io::IO, ::MIME"text/plain", result::ReducedFrequencyCollision)
    _show_stage_plain(io, "Reduced causal-frequency collision", result)
    write(io, "\n  I_reg(k) = ")
    _write_plain_statistical_collision(io, reduced_regular_terms(result))
    print(
        io,
        "\n  causal blockers: ",
        length(reduced_blocked_terms(result)),
        "\n  unresolved Trotter states: ",
        length(reduced_trotter_terms(result)),
    )
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", result::ReducedFrequencyCollision)
    return _latex_display(io) do
        write(io, "\\begin{aligned}\nI_{\\mathrm{reg}}(k)={}&")
        _write_latex_statistical_collision(io, reduced_regular_terms(result))
        write(
            io,
            "\\\\\nN_{\\mathrm{blocked}}={}&",
            string(length(reduced_blocked_terms(result))),
            ",\\qquad N_{\\mathrm{Trotter}}=",
            string(length(reduced_trotter_terms(result))),
            "\n\\end{aligned}",
        )
        return nothing
    end
end

function Base.show(io::IO, ::MIME"text/plain", result::OccupationReducedExpression)
    _show_stage_plain(io, "Occupation-reduced regular collision", result)
    write(io, "\n  C_n^reg(k) = ")
    _write_plain_occupation_collision(io, occupation_reduced_terms(result))
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", result::OccupationReducedExpression)
    return _latex_display(io) do
        write(io, "\\begin{aligned}\nC_n^{\\mathrm{reg}}(k)={}&")
        _write_latex_occupation_collision(io, occupation_reduced_terms(result))
        write(io, "\n\\end{aligned}")
        return nothing
    end
end

function Base.show(io::IO, ::MIME"text/plain", result::LoopQuotientedExpression)
    _show_stage_plain(io, "Loop-quotiented regular collision", result)
    write(io, "\n  C_n^quot(k) = ")
    _write_plain_occupation_collision(io, loop_quotient_terms(result))
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", result::LoopQuotientedExpression)
    return _latex_display(io) do
        write(io, "\\begin{aligned}\nC_n^{\\mathrm{quot}}(k)={}&")
        _write_latex_occupation_collision(io, loop_quotient_terms(result))
        write(io, "\n\\end{aligned}")
        return nothing
    end
end

function Base.show(io::IO, ::MIME"text/plain", kernel::CollisionKernel)
    _show_stage_plain(io, "Collision kernel", kernel)
    write(io, "\n  C_n(k) = ")
    _write_plain_occupation_collision(io, collision_kernel_terms(kernel))
    return nothing
end
function Base.show(io::IO, ::MIME"text/latex", kernel::CollisionKernel)
    return _latex_display(io) do
        write(io, "\\begin{aligned}\nC_n(k)={}&")
        _write_latex_occupation_collision(io, collision_kernel_terms(kernel))
        write(io, "\n\\end{aligned}")
        return nothing
    end
end
