# Exact display-only compaction for physical collision expressions.
#
# This deliberately does not mutate compiler IR. The collision IR is first regrouped by
# physical sector, then the kinematic/statistical tensor product is reconstructed exactly.
# Small low-degree kinematic polynomials are matched against products of routed linear forms,
# while statistical and occupation polynomials are transformed to the physically natural
# tensor-product bases (1±F) and n/(1+σn). The shorter exact representation is rendered.

function _inline_sum_string(terms::Vector{String})
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

function _wrap_sum(text::String, nterms::Int, latex::Bool)
    nterms <= 1 && return text
    return latex ? "\\left[" * text * "\\right]" : "[" * text * "]"
end

_display_cost(text::String) = ncodeunits(text)

function _expanded_polynomial_string(polynomial, basis, external, latex::Bool)
    rendered = String[]
    for (monomial, coefficient) in polynomial
        monomial_text = _polynomial_monomial_string(monomial, basis, external, latex)
        factors = isempty(monomial_text) ? String[] : String[monomial_text]
        push!(rendered, _term_string(coefficient, factors, latex))
    end
    return _wrap_sum(_inline_sum_string(rendered), length(rendered), latex)
end

function _collect_polynomial_atoms(polynomial)
    atoms = Any[]
    for (monomial, _) in polynomial
        for atom in monomial
            any(existing -> isequal(existing, atom), atoms) || push!(atoms, atom)
        end
    end
    sort!(atoms)
    return atoms
end

function _multilinear_coefficients(polynomial, atoms)
    n = length(atoms)
    n <= 8 || return nothing
    atom_index = Dict{Any,Int}(atom => i for (i, atom) in enumerate(atoms))
    C = isempty(polynomial) ? Rational{Int} : typeof(last(first(polynomial)))
    coefficients = fill(zero(C), 1 << n)
    for (monomial, coefficient) in polynomial
        mask = 0
        for atom in monomial
            index = atom_index[atom]
            bit = 1 << (index - 1)
            iszero(mask & bit) || return nothing
            mask |= bit
        end
        coefficients[mask + 1] += coefficient
    end
    return coefficients
end

function _occupation_literal_string(
    atom, gain::Bool, basis, external, latex::Bool, sigma::Int
)
    atom_text = _distribution_atom_string(
        "n", atom.family, atom.momentum, basis, external, latex
    )
    gain || return atom_text
    sign = sigma > 0 ? "+" : "-"
    return if latex
        "\\left(1" * sign * atom_text * "\\right)"
    else
        "(1" * sign * atom_text * ")"
    end
end

function _statistical_literal_string(atom, plus::Bool, basis, external, latex::Bool)
    atom_text = _distribution_atom_string(
        "F", atom.family, atom.momentum, basis, external, latex
    )
    sign = plus ? "+" : "-"
    return if latex
        "\\left(1" * sign * atom_text * "\\right)"
    else
        "(1" * sign * atom_text * ")"
    end
end

function _physical_basis_occupation_string(
    polynomial::OccupationPolynomial{C,S}, basis, external, latex::Bool
) where {C<:Number,S<:Statistics}
    atoms = _collect_polynomial_atoms(polynomial)
    coefficients = _multilinear_coefficients(polynomial, atoms)
    coefficients === nothing && return nothing
    n = length(atoms)
    sigma = Int(occupation_statistics_sign(S))

    # Per variable: a0 + a1 n = c_n n + c_g (1 + σ n),
    # c_n = a1 - σ a0, c_g = a0.
    for variable in 0:(n - 1)
        bit = 1 << variable
        for mask in 0:((1 << n) - 1)
            iszero(mask & bit) || continue
            i0 = mask + 1
            i1 = (mask | bit) + 1
            a0 = coefficients[i0]
            a1 = coefficients[i1]
            coefficients[i0] = a1 - sigma * a0
            coefficients[i1] = a0
        end
    end

    rendered = String[]
    for state in 0:((1 << n) - 1)
        coefficient = coefficients[state + 1]
        iszero(coefficient) && continue
        factors = String[]
        for (index, atom) in enumerate(atoms)
            gain = !iszero(state & (1 << (index - 1)))
            push!(
                factors,
                _occupation_literal_string(atom, gain, basis, external, latex, sigma),
            )
        end
        push!(rendered, _term_string(coefficient, factors, latex))
    end
    sort!(rendered; by=text -> (startswith(text, "-"), text))
    return _wrap_sum(_inline_sum_string(rendered), length(rendered), latex)
end

function _physical_basis_statistical_string(
    polynomial::StatisticalPolynomial{C,S}, basis, external, latex::Bool
) where {C<:Number,S<:Statistics}
    atoms = _collect_polynomial_atoms(polynomial)
    raw = _multilinear_coefficients(polynomial, atoms)
    raw === nothing && return nothing
    n = length(atoms)
    D = promote_type(C, Rational{Int})
    coefficients = D[convert(D, value) for value in raw]

    # Per variable: a0 + a1 F = c_- (1-F) + c_+ (1+F),
    # c_-=(a0-a1)/2, c_+=(a0+a1)/2.
    for variable in 0:(n - 1)
        bit = 1 << variable
        for mask in 0:((1 << n) - 1)
            iszero(mask & bit) || continue
            i0 = mask + 1
            i1 = (mask | bit) + 1
            a0 = coefficients[i0]
            a1 = coefficients[i1]
            coefficients[i0] = (a0 - a1) / 2
            coefficients[i1] = (a0 + a1) / 2
        end
    end

    rendered = String[]
    for state in 0:((1 << n) - 1)
        coefficient = coefficients[state + 1]
        iszero(coefficient) && continue
        factors = String[]
        for (index, atom) in enumerate(atoms)
            plus = !iszero(state & (1 << (index - 1)))
            push!(factors, _statistical_literal_string(atom, plus, basis, external, latex))
        end
        push!(rendered, _term_string(coefficient, factors, latex))
    end
    sort!(rendered; by=text -> (startswith(text, "-"), text))
    return _wrap_sum(_inline_sum_string(rendered), length(rendered), latex)
end

function _compact_distribution_string(polynomial, basis, external, latex::Bool)
    expanded = _expanded_polynomial_string(polynomial, basis, external, latex)
    physical = if polynomial isa OccupationPolynomial
        _physical_basis_occupation_string(polynomial, basis, external, latex)
    elseif polynomial isa StatisticalPolynomial
        _physical_basis_statistical_string(polynomial, basis, external, latex)
    else
        nothing
    end
    physical === nothing && return expanded
    return _display_cost(physical) < _display_cost(expanded) ? physical : expanded
end

function _real_rational(value)
    value isa Rational && return value
    if value isa Complex && iszero(imag(value)) && real(value) isa Rational
        return real(value)
    end
    return nothing
end

function _polynomial_content(polynomial)
    isempty(polynomial) && return 1 // 1
    rationals = Rational{Int}[]
    for (_, coefficient) in polynomial
        value = _real_rational(coefficient)
        value === nothing && return 1 // 1
        push!(rationals, convert(Rational{Int}, value))
    end
    numerator_gcd = zero(Int)
    denominator_lcm = one(Int)
    for value in rationals
        numerator_gcd = gcd(numerator_gcd, abs(numerator(value)))
        denominator_lcm = lcm(denominator_lcm, denominator(value))
    end
    iszero(numerator_gcd) && return 1 // 1
    content = numerator_gcd // denominator_lcm
    return first(rationals) < 0 ? -content : content
end

function _primitive_polynomial(polynomial)
    content = _polynomial_content(polynomial)
    isone(content) && return content, polynomial
    return content, inv(content) * polynomial
end

function _polynomial_scale(a, b)
    length(a) == length(b) || return nothing
    isempty(a) && return one(Rational{Int})
    left = collect(a)
    right = collect(b)
    for i in eachindex(left)
        isequal(first(left[i]), first(right[i])) || return nothing
    end
    scale = last(first(left)) / last(first(right))
    for i in eachindex(left)
        last(left[i]) == scale * last(right[i]) || return nothing
    end
    return scale
end

function _basis_component_polynomial(component::MomentumComponent, basis::MomentumBasis)
    terms = Pair{MomentumMonomial,ComplexRationals}[]
    for slot in eachindex(component.momentum.coefficients)
        coefficient = component.momentum[slot]
        iszero(coefficient) && continue
        unit = MomentumComponent(basis_momentum(basis, slot), component.axis)
        monomial = MomentumMonomial(MomentumComponent[unit])
        push!(terms, monomial => convert(ComplexRationals, coefficient))
    end
    return MomentumPolynomial{ComplexRationals}(terms)
end

function _expanded_momentum_monomial(monomial::MomentumMonomial, basis::MomentumBasis)
    out = MomentumPolynomial(MomentumMonomial(), one(ComplexRationals))
    for component in monomial
        out *= _basis_component_polynomial(component, basis)
    end
    return out
end

function _expanded_kinematic_polynomial(terms, basis::MomentumBasis)
    out = MomentumPolynomial{ComplexRationals}()
    for (monomial, coefficient) in terms
        converted = try
            convert(ComplexRationals, coefficient)
        catch
            return nothing
        end
        out += converted * _expanded_momentum_monomial(monomial, basis)
    end
    return out
end

function _normalize_candidate_momentum(momentum_value::LinearMomentum)
    pivot = findfirst(value -> !iszero(value), momentum_value.coefficients)
    pivot === nothing && return nothing
    return inv(momentum_value[pivot]) * momentum_value
end

function _push_unique_momentum!(
    momenta::Vector{LinearMomentum}, momentum_value::LinearMomentum
)
    normalized = _normalize_candidate_momentum(momentum_value)
    normalized === nothing && return momenta
    any(value -> isequal(value, normalized), momenta) || push!(momenta, normalized)
    return momenta
end

function _candidate_momenta(rows, basis::MomentumBasis)
    raw = LinearMomentum[]
    for slot in 1:length(basis)
        push!(raw, basis_momentum(basis, slot))
    end
    for monomial in sort!(collect(keys(rows)))
        for component in monomial
            push!(raw, component.momentum)
        end
    end

    candidates = LinearMomentum[]
    for momentum_value in raw
        _push_unique_momentum!(candidates, momentum_value)
    end
    nraw = length(raw)
    for i in 1:nraw, j in (i + 1):nraw
        j <= nraw || continue
        _push_unique_momentum!(candidates, raw[i] - raw[j])
    end
    return candidates
end

function _axis_counts(polynomial::MomentumPolynomial)
    isempty(polynomial) && return Dict{Symbol,Int}()
    reference = Dict{Symbol,Int}()
    first_monomial = first(first(polynomial))
    for component in first_monomial
        reference[component.axis] = get(reference, component.axis, 0) + 1
    end
    for (monomial, _) in polynomial
        counts = Dict{Symbol,Int}()
        for component in monomial
            counts[component.axis] = get(counts, component.axis, 0) + 1
        end
        counts == reference || return nothing
    end
    return reference
end

function _momentum_polynomial_scale(
    target::MomentumPolynomial, candidate::MomentumPolynomial
)
    length(target) == length(candidate) || return nothing
    isempty(target) && return one(ComplexRationals)
    target_terms = collect(target)
    candidate_terms = collect(candidate)
    for i in eachindex(target_terms)
        isequal(first(target_terms[i]), first(candidate_terms[i])) || return nothing
    end
    scale = last(first(target_terms)) / last(first(candidate_terms))
    for i in eachindex(target_terms)
        last(target_terms[i]) == scale * last(candidate_terms[i]) || return nothing
    end
    return scale
end

function _factor_components_string(factors, basis, external, latex::Bool)
    rendered = _grouped_factor_strings(factors) do component
        return _component_string(component, basis, external, latex)
    end
    return join(rendered, latex ? "\\," : " ")
end

function _search_linear_factorization(
    target::MomentumPolynomial,
    candidates::Vector{MomentumComponent},
    counts::Dict{Symbol,Int},
    basis,
    external,
    latex::Bool,
)
    degree = sum(values(counts))
    degree <= 6 || return nothing
    max_candidates = degree <= 4 ? 24 : 14
    length(candidates) > max_candidates && resize!(candidates, max_candidates)

    remaining = copy(counts)
    identity_poly = MomentumPolynomial(MomentumMonomial(), one(ComplexRationals))
    factors = MomentumComponent[]
    best = Ref{Union{Nothing,Tuple{ComplexRationals,String}}}(nothing)

    function recurse(start::Int, depth::Int, product::MomentumPolynomial{ComplexRationals})
        if depth > degree
            scale = _momentum_polynomial_scale(target, product)
            scale === nothing && return nothing
            text = _factor_components_string(factors, basis, external, latex)
            candidate = (scale, text)
            if best[] === nothing || _display_cost(text) < _display_cost(last(best[]))
                best[] = candidate
            end
            return nothing
        end
        for index in start:length(candidates)
            component = candidates[index]
            available = get(remaining, component.axis, 0)
            iszero(available) && continue
            remaining[component.axis] = available - 1
            push!(factors, component)
            recurse(
                index, depth + 1, product * _basis_component_polynomial(component, basis)
            )
            pop!(factors)
            remaining[component.axis] = available
        end
        return nothing
    end

    recurse(1, 1, identity_poly)
    return best[]
end

function _expanded_kinematic_string(terms, basis, external, latex::Bool)
    rendered = String[]
    for (monomial, coefficient) in terms
        monomial_text = _momentum_monomial_string(monomial, basis, external, latex)
        factors = isempty(monomial_text) ? String[] : String[monomial_text]
        push!(rendered, _term_string(coefficient, factors, latex))
    end
    return _wrap_sum(_inline_sum_string(rendered), length(rendered), latex)
end

function _compact_kinematic_string(terms, rows, basis, external, latex::Bool)
    fallback_content = one(Rational{Int})
    fallback_terms = collect(terms)
    if !isempty(fallback_terms)
        rational_candidates = [_real_rational(last(term)) for term in fallback_terms]
        if all(value -> value !== nothing, rational_candidates)
            rationals = Rational{Int}[
                convert(Rational{Int}, value) for value in rational_candidates
            ]
            numerator_gcd = foldl(
                gcd, (abs(numerator(value)) for value in rationals); init=0
            )
            denominator_lcm = foldl(
                lcm, (denominator(value) for value in rationals); init=1
            )
            if !iszero(numerator_gcd)
                fallback_content = numerator_gcd // denominator_lcm
                first(rationals) < 0 && (fallback_content = -fallback_content)
                fallback_terms = [
                    first(term) => last(term) / fallback_content for term in fallback_terms
                ]
            end
        end
    end
    fallback = _expanded_kinematic_string(fallback_terms, basis, external, latex)

    expanded = _expanded_kinematic_polynomial(terms, basis)
    expanded === nothing && return fallback_content, fallback
    counts = _axis_counts(expanded)
    counts === nothing && return fallback_content, fallback
    degree = sum(values(counts))
    iszero(degree) && return fallback_content, ""

    momenta = _candidate_momenta(rows, basis)
    axes = sort!(collect(keys(counts)))
    candidates = MomentumComponent[]
    for momentum_value in momenta, axis in axes
        push!(candidates, MomentumComponent(momentum_value, axis))
    end
    factorization = _search_linear_factorization(
        expanded, candidates, counts, basis, external, latex
    )
    factorization === nothing && return fallback_content, fallback

    factor_scale, factor_text = factorization
    factor_render = _term_string(factor_scale, String[factor_text], latex)
    fallback_render = _term_string(fallback_content, String[fallback], latex)
    if _display_cost(factor_render) < _display_cost(fallback_render)
        return factor_scale, factor_text
    end
    return fallback_content, fallback
end

function _collect_kinematic_rows(group)
    rows = Dict{MomentumMonomial,Any}()
    for (sector, polynomial) in group
        for (kinematic, kinematic_coefficient) in kinematic_factor(sector)
            contribution = kinematic_coefficient * polynomial
            if haskey(rows, kinematic)
                combined = rows[kinematic] + contribution
                if iszero(combined)
                    delete!(rows, kinematic)
                else
                    rows[kinematic] = combined
                end
            elseif !iszero(contribution)
                rows[kinematic] = contribution
            end
        end
    end
    return rows
end

function _rank_one_rows(rows)
    isempty(rows) && return nothing
    row_keys = sort!(collect(keys(rows)))
    reference = rows[first(row_keys)]
    scales = Dict{MomentumMonomial,Any}()
    for key in row_keys
        scale = _polynomial_scale(rows[key], reference)
        scale === nothing && return nothing
        scales[key] = scale
    end
    content, primitive = _primitive_polynomial(reference)
    kinematic = Pair{MomentumMonomial,Any}[key => scales[key] * content for key in row_keys]
    return kinematic, primitive
end

function _compact_group_string(group, latex::Bool)
    sector = first(first(group))
    basis = momentum_basis(sector)
    external = external_wigner_momentum(sector)
    parameter_text = _parameter_string(parameters(sector), latex)
    measure_text = _measure_string(basis, external, latex)
    support_text = _frequency_support_string(
        frequency_support(sector), basis, external, latex
    )
    rows = _collect_kinematic_rows(group)
    factorization = _rank_one_rows(rows)

    if factorization !== nothing
        kinematic_terms, distribution = factorization
        coefficient, kinematic_text = _compact_kinematic_string(
            kinematic_terms, rows, basis, external, latex
        )
        distribution_text = _compact_distribution_string(
            distribution, basis, external, latex
        )
        factors = String[
            parameter_text,
            measure_text,
            kinematic_text,
            distribution_text == "1" ? "" : distribution_text,
            support_text,
        ]
        return _term_string(coefficient, factors, latex)
    end

    rendered = String[]
    for key in sort!(collect(keys(rows)))
        content, primitive = _primitive_polynomial(rows[key])
        coefficient, kinematic_text = _compact_kinematic_string(
            [key => content], rows, basis, external, latex
        )
        distribution_text = _compact_distribution_string(primitive, basis, external, latex)
        factors = String[
            parameter_text,
            measure_text,
            kinematic_text,
            distribution_text == "1" ? "" : distribution_text,
            support_text,
        ]
        push!(rendered, _term_string(coefficient, factors, latex))
    end
    return _inline_sum_string(rendered)
end

function _compact_physical_collision_string(terms, latex::Bool)
    groups = Dict{Any,Vector{Any}}()
    for (sector, polynomial) in terms
        key = (
            parameters(sector),
            momentum_basis(sector),
            external_wigner_momentum(sector),
            frequency_support(sector),
        )
        push!(get!(groups, key, Any[]), sector => polynomial)
    end
    rendered = String[]
    for group in values(groups)
        push!(rendered, _compact_group_string(group, latex))
    end
    sort!(rendered)
    return _sum_string(rendered, latex)
end

function Base.show(
    io::IO, ::MIME"text/plain", result::ReducedFrequencyCollision{C,S,O,E1,E2,G,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    _show_stage_plain(io, "Reduced causal-frequency collision", result)
    write(
        io,
        "\n  I_reg(k) = ",
        _compact_physical_collision_string(reduced_regular_terms(result), false),
    )
    write(io, "\n  ", _blocked_summary(result, false))
    return nothing
end

function Base.show(
    io::IO, ::MIME"text/latex", result::ReducedFrequencyCollision{C,S,O,E1,E2,G,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    return _latex_display(io) do
        write(
            io,
            "\\begin{aligned}\nI_{\\mathrm{reg}}(k)={}&",
            _compact_physical_collision_string(reduced_regular_terms(result), true),
            "\\\\\n&",
            _blocked_summary(result, true),
            "\n\\end{aligned}",
        )
        return nothing
    end
end

function Base.show(
    io::IO, ::MIME"text/plain", result::OccupationReducedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,G,Ctx<:AbstractWignerContext}
    _show_stage_plain(io, "Occupation-reduced regular collision", result)
    write(
        io,
        "\n  C_n^reg(k) = ",
        _compact_physical_collision_string(occupation_reduced_terms(result), false),
    )
    return nothing
end

function Base.show(
    io::IO, ::MIME"text/latex", result::OccupationReducedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,G,Ctx<:AbstractWignerContext}
    return _latex_display(io) do
        write(
            io,
            "\\begin{aligned}\nC_n^{\\mathrm{reg}}(k)={}&",
            _compact_physical_collision_string(occupation_reduced_terms(result), true),
            "\n\\end{aligned}",
        )
        return nothing
    end
end

function Base.show(
    io::IO, ::MIME"text/plain", result::LoopQuotientedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,G,Ctx<:AbstractWignerContext}
    _show_stage_plain(io, "Loop-quotiented regular collision", result)
    write(
        io,
        "\n  C_n^quot(k) = ",
        _compact_physical_collision_string(loop_quotient_terms(result), false),
    )
    return nothing
end

function Base.show(
    io::IO, ::MIME"text/latex", result::LoopQuotientedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,G,Ctx<:AbstractWignerContext}
    return _latex_display(io) do
        write(
            io,
            "\\begin{aligned}\nC_n^{\\mathrm{quot}}(k)={}&",
            _compact_physical_collision_string(loop_quotient_terms(result), true),
            "\n\\end{aligned}",
        )
        return nothing
    end
end

function Base.show(
    io::IO, ::MIME"text/plain", kernel::CollisionKernel{C,S,O,G,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,G,Ctx<:AbstractWignerContext}
    _show_stage_plain(io, "Collision kernel", kernel)
    write(
        io,
        "\n  C_n(k) = ",
        _compact_physical_collision_string(collision_kernel_terms(kernel), false),
    )
    return nothing
end

function Base.show(
    io::IO, ::MIME"text/latex", kernel::CollisionKernel{C,S,O,G,Ctx}
) where {C<:Number,S<:_PhysicalDisplayStatistics,O,G,Ctx<:AbstractWignerContext}
    return _latex_display(io) do
        write(
            io,
            "\\begin{aligned}\nC_n(k)={}&",
            _compact_physical_collision_string(collision_kernel_terms(kernel), true),
            "\n\\end{aligned}",
        )
        return nothing
    end
end