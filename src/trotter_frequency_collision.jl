"""
Frequency-independent identity of a post-Trotter causal expression.

`sector` is the same topology-free physical identity used by `CanonicalFrequencyCollision`.
`active_frequency_indices` records the original canonical loop-frequency coordinates that still
require contour integration. Contributions are combined only when both pieces agree, so terms
are never compared after silently changing frequency stages.
"""
struct TrotterFrequencyGroupKey{S<:Statistics}
    sector::CanonicalFrequencySector{S}
    active_frequency_indices::Vector{Int}
end

function TrotterFrequencyGroupKey(
    sector::CanonicalFrequencySector{S}, active_frequency_indices::AbstractVector{Int}
) where {S<:Statistics}
    active = collect(active_frequency_indices)
    sort!(active)
    unique!(active)
    return TrotterFrequencyGroupKey{S}(sector, active)
end

function Base.isequal(
    a::TrotterFrequencyGroupKey{S}, b::TrotterFrequencyGroupKey{S}
) where {S<:Statistics}
    return isequal(a.sector, b.sector) &&
           isequal(a.active_frequency_indices, b.active_frequency_indices)
end
Base.:(==)(a::TrotterFrequencyGroupKey, b::TrotterFrequencyGroupKey) = isequal(a, b)
function Base.hash(key::TrotterFrequencyGroupKey, h::UInt)
    return hash(
        Tuple(key.active_frequency_indices),
        hash(key.sector, hash(TrotterFrequencyGroupKey, h)),
    )
end

canonical_frequency_sector(key::TrotterFrequencyGroupKey) = key.sector
function active_trotter_frequencies(key::TrotterFrequencyGroupKey)
    return copy(key.active_frequency_indices)
end

"""
Canonical result of all locally certified equal-time integrations in the shifted sector.

`expressions` are complete causal sums grouped before any contour decision. `constants` contain
sectors for which Trotter integration removed every loop frequency exactly. `unresolved` retains
partially reduced states whose remaining shifted factors are not structurally local and hence
cannot be assigned a universal equal-time value.
"""
struct TrotterFrequencyCollision{C<:Number,S<:Statistics,E1,E2}
    expressions::Dict{TrotterFrequencyGroupKey{S},CausalFrequencyExpression{C,S}}
    constants::Dict{CanonicalFrequencySector{S},C}
    unresolved::Vector{TrotterFrequencyState{C,S,E1,E2}}
end

function grouped_trotter_expressions(collision::TrotterFrequencyCollision)
    return collision.expressions
end
function grouped_trotter_constants(collision::TrotterFrequencyCollision)
    return collision.constants
end
function unresolved_trotter_states(collision::TrotterFrequencyCollision)
    return collision.unresolved
end
function Base.isempty(collision::TrotterFrequencyCollision)
    return isempty(collision.expressions) &&
           isempty(collision.constants) &&
           isempty(collision.unresolved)
end

"""One canonical causal blocker after ordinary and post-Trotter stages have been merged."""
struct TrotterFrequencyBlockedContribution{C<:Number,S<:Statistics}
    key::TrotterFrequencyGroupKey{S}
    expression::CausalFrequencyExpression{C,S}
    kind::CausalFrequencyIntegrationKind
    blocked_frequency::Int
    supports::Vector{AffineSingularSupport{S}}
end

function canonical_frequency_sector(blocked::TrotterFrequencyBlockedContribution)
    return canonical_frequency_sector(blocked.key)
end
function active_trotter_frequencies(blocked::TrotterFrequencyBlockedContribution)
    return active_trotter_frequencies(blocked.key)
end
function trotter_frequency_blocker_kind(blocked::TrotterFrequencyBlockedContribution)
    return blocked.kind
end
function trotter_frequency_blocked_index(blocked::TrotterFrequencyBlockedContribution)
    return blocked.blocked_frequency
end
function trotter_frequency_blocked_expression(blocked::TrotterFrequencyBlockedContribution)
    return blocked.expression
end
function trotter_frequency_blocked_supports(blocked::TrotterFrequencyBlockedContribution)
    return blocked.supports
end

"""
Final #303 reduction state after exact equal-time lowering and staged canonical contour reduction.

Every ordinary and shifted contribution of the same physical sector is brought to a common
active-frequency stage before the next contour decision. `boundary` retains canonical expressions
with no remaining loop-frequency integration but nontrivial causal energy denominators. These are
terminal #303 objects and the direct input of the later repeated-pole/Plemelj boundary-value layer.
`blocked` contains genuine frequency-integration blockers of the complete expression available at
that stage. Boundary-value evaluation itself remains deferred.
"""
struct CanonicalTrotterFrequencyReduction{C<:Number,S<:Statistics,E1,E2}
    constants::Dict{CanonicalFrequencySector{S},C}
    boundary::Dict{CanonicalFrequencySector{S},CausalFrequencyExpression{C,S}}
    blocked::Dict{TrotterFrequencyGroupKey{S},TrotterFrequencyBlockedContribution{C,S}}
    unresolved::Vector{TrotterFrequencyState{C,S,E1,E2}}
end

function canonical_trotter_constants(reduction::CanonicalTrotterFrequencyReduction)
    return reduction.constants
end
function canonical_trotter_boundary_expressions(
    reduction::CanonicalTrotterFrequencyReduction
)
    return reduction.boundary
end
function blocked_trotter_contributions(reduction::CanonicalTrotterFrequencyReduction)
    return reduction.blocked
end
function unresolved_trotter_states(reduction::CanonicalTrotterFrequencyReduction)
    return reduction.unresolved
end
function Base.isempty(reduction::CanonicalTrotterFrequencyReduction)
    return isempty(reduction.constants) &&
           isempty(reduction.boundary) &&
           isempty(reduction.blocked) &&
           isempty(reduction.unresolved)
end

function _complex_trotter_frequency_state(
    state::TrotterFrequencyState{C,S,E1,E2}
) where {C<:Number,S<:Statistics,E1,E2}
    D = promote_type(C, ComplexRationals)
    return TrotterFrequencyState{D,S,E1,E2}(
        state.source,
        convert(D, state.coefficient),
        state.statistical,
        state.parameter,
        copy(state.active_line_indices),
        copy(state.active_frequency_indices),
    )
end

function _trotter_frequency_sector(
    state::TrotterFrequencyState{C,S}
) where {C,S<:Statistics}
    return _canonical_frequency_sector(state.source, state.statistical, state.parameter)
end

function _merge_trotter_causal_expression(
    a::CausalFrequencyExpression{C,S}, b::CausalFrequencyExpression{C,S}
) where {C<:Number,S<:Statistics}
    terms = CausalFrequencyTerm{C,S}[]
    sizehint!(terms, length(a.terms) + length(b.terms))
    append!(terms, a.terms)
    append!(terms, b.terms)
    return CausalFrequencyExpression(terms)
end

function _merge_trotter_stage!(
    stages::Dict{TrotterFrequencyGroupKey{S},CausalFrequencyExpression{C,S}},
    key::TrotterFrequencyGroupKey{S},
    expression::CausalFrequencyExpression{C,S},
) where {C<:Number,S<:Statistics}
    isempty(expression) && return stages
    if haskey(stages, key)
        merged = _merge_trotter_causal_expression(stages[key], expression)
        if isempty(merged)
            delete!(stages, key)
        else
            stages[key] = merged
        end
    else
        stages[key] = expression
    end
    return stages
end

function _constant_trotter_expression(
    coefficient::C, ::CanonicalFrequencySector{S}
) where {C<:Number,S<:Statistics}
    return CausalFrequencyExpression(
        CausalFrequencyTerm{C,S}[CausalFrequencyTerm(
            coefficient, CausalFrequencyDenominator{S}[]
        )],
    )
end

function _causal_frequency_dimension(expression::CausalFrequencyExpression)
    dimension = 0
    for term in causal_frequency_terms(expression)
        for denominator in causal_frequency_denominators(term)
            candidate = length(denominator.loop_coefficients)
            if iszero(dimension)
                dimension = candidate
            elseif candidate != dimension
                throw(
                    DimensionMismatch(
                        "causal denominators use different frequency dimensions"
                    ),
                )
            end
        end
    end
    return dimension
end

function _same_trotter_sector(a::TrotterFrequencyGroupKey, b::TrotterFrequencyGroupKey)
    return isequal(canonical_frequency_sector(a), canonical_frequency_sector(b))
end

function _trotter_stage_subset(
    candidate::TrotterFrequencyGroupKey, source::TrotterFrequencyGroupKey
)
    _same_trotter_sector(candidate, source) || return false
    length(candidate.active_frequency_indices) < length(source.active_frequency_indices) ||
        return false
    return all(
        index -> index in source.active_frequency_indices,
        candidate.active_frequency_indices,
    )
end

function _trotter_active_indices_less(a::Vector{Int}, b::Vector{Int})
    for index in eachindex(a, b)
        a[index] == b[index] || return a[index] < b[index]
    end
    return length(a) < length(b)
end

function _nearest_existing_trotter_subset(
    stages::Dict{TrotterFrequencyGroupKey{S},CausalFrequencyExpression{C,S}},
    source::TrotterFrequencyGroupKey{S},
) where {C<:Number,S<:Statistics}
    best = source
    best_length = -1
    found = false
    for candidate in keys(stages)
        _trotter_stage_subset(candidate, source) || continue
        candidate_length = length(candidate.active_frequency_indices)
        if !found ||
            candidate_length > best_length ||
            (
                candidate_length == best_length && _trotter_active_indices_less(
                    candidate.active_frequency_indices, best.active_frequency_indices
                )
            )
            best = candidate
            best_length = candidate_length
            found = true
        end
    end
    return found, best
end

function _trotter_sector_common_frequencies(
    stages::Dict{TrotterFrequencyGroupKey{S},CausalFrequencyExpression{C,S}},
    source::TrotterFrequencyGroupKey{S},
) where {C<:Number,S<:Statistics}
    common = copy(source.active_frequency_indices)
    found_other = false
    for candidate in keys(stages)
        candidate === source && continue
        _same_trotter_sector(candidate, source) || continue
        found_other = true
        filter!(index -> index in candidate.active_frequency_indices, common)
    end
    return found_other ? common : Int[]
end

function _next_trotter_frequency(
    stages::Dict{TrotterFrequencyGroupKey{S},CausalFrequencyExpression{C,S}},
    source::TrotterFrequencyGroupKey{S},
) where {C<:Number,S<:Statistics}
    active = source.active_frequency_indices
    isempty(active) && return 0

    found_subset, subset = _nearest_existing_trotter_subset(stages, source)
    if found_subset
        for index in active
            index in subset.active_frequency_indices || return index
        end
    end

    common = _trotter_sector_common_frequencies(stages, source)
    for index in active
        index in common || return index
    end
    return first(active)
end

function _maximal_trotter_stage(
    stages::Dict{TrotterFrequencyGroupKey{S},CausalFrequencyExpression{C,S}}
) where {C<:Number,S<:Statistics}
    isempty(stages) &&
        throw(ArgumentError("cannot select a stage from an empty dictionary"))
    selected = first(keys(stages))
    selected_length = length(selected.active_frequency_indices)
    for key in keys(stages)
        candidate_length = length(key.active_frequency_indices)
        candidate_length > selected_length || continue
        selected = key
        selected_length = candidate_length
    end
    return selected
end

function _is_trotter_constant_expression(expression::CausalFrequencyExpression)
    terms = causal_frequency_terms(expression)
    length(terms) == 1 || return false
    return isempty(causal_frequency_denominators(only(terms)))
end

function _trotter_constant_coefficient(
    expression::CausalFrequencyExpression{C,S}
) where {C,S}
    _is_trotter_constant_expression(expression) ||
        throw(ArgumentError("causal expression is not a scalar constant"))
    return causal_frequency_coefficient(only(causal_frequency_terms(expression)))
end

"""
    reduce_shifted_trotter_frequencies(collision)

Apply every structurally certified local equal-time rule to the shifted part of one canonical
collision and then recombine the residual algebra before contour reduction.

This operation is deliberately not a second frequency integrator. It only removes isolated
finite-Trotter frequencies, preserves the original spatial momentum basis, and returns residual
`CausalFrequencyExpression`s in the same representation consumed by #299/#301. No topology,
perturbative order, or graph label enters the grouping key.
"""
function reduce_shifted_trotter_frequencies(
    collision::CanonicalFrequencyCollision{C,S,O,E1,E2}
) where {C<:Number,S<:Statistics,O,E1,E2}
    D = promote_type(C, ComplexRationals)
    expressions = Dict{TrotterFrequencyGroupKey{S},CausalFrequencyExpression{D,S}}()
    constants = Dict{CanonicalFrequencySector{S},D}()
    unresolved = TrotterFrequencyState{D,S,E1,E2}[]

    for term in shifted_frequency_terms(collision)
        state = _complex_trotter_frequency_state(TrotterFrequencyState(term))
        reduced = eliminate_isolated_trotter_frequencies(state)
        if has_unresolved_trotter_frequency(reduced)
            push!(unresolved, reduced)
            continue
        end

        sector = _trotter_frequency_sector(reduced)
        if trotter_frequency_complete(reduced)
            constants[sector] =
                get(constants, sector, zero(D)) + convert(D, source_coefficient(reduced))
            iszero(constants[sector]) && delete!(constants, sector)
            continue
        end

        key = TrotterFrequencyGroupKey(sector, active_trotter_frequencies(reduced))
        residual = trotter_residual_causal_expression(reduced, target_family(collision))
        _merge_trotter_stage!(expressions, key, residual)
    end

    return TrotterFrequencyCollision{D,S,E1,E2}(expressions, constants, unresolved)
end

"""
    reduce_canonical_trotter_frequencies(collision)

Compile the complete canonical collision through the exact equal-time layer without splitting
ordinary and Trotter-reduced contributions at a later common frequency stage.

Ordinary causal sectors enter at their full loop-frequency set. Shifted terms first undergo only
certified local equal-time elimination and enter at the active set that remains. The reducer then
advances the highest active stage by one canonical frequency at a time through frozen #299/#301
semantics. Whenever another contribution already exists at the resulting stage, exact canonical
recombination occurs before any further contour decision. This is the collision-level analogue
of SSA/data-flow joining: a contour classification is made only after every contribution that
can reach the same frequency stage has been merged.

Expressions with no remaining loop-frequency integration are retained canonically. Pure scalars
are separated into `constants`; nontrivial causal energy-boundary expressions are preserved for
#306. Deferred zero-prescription, repeated-pole and pinch results that still block a loop-frequency
integration remain typed blockers. No finite regulator width or topology-specific rule is used.
"""
function reduce_canonical_trotter_frequencies(
    collision::CanonicalFrequencyCollision{C,S,O,E1,E2}
) where {C<:Number,S<:Statistics,O,E1,E2}
    D = promote_type(C, ComplexRationals)
    shifted = reduce_shifted_trotter_frequencies(collision)
    stages = Dict{TrotterFrequencyGroupKey{S},CausalFrequencyExpression{D,S}}()

    for (sector, expression) in canonical_frequency_expressions(collision)
        canonical = _complex_causal_frequency_expression(expression)
        dimension = _causal_frequency_dimension(canonical)
        key = TrotterFrequencyGroupKey(sector, collect(1:dimension))
        _merge_trotter_stage!(stages, key, canonical)
    end
    for (key, expression) in grouped_trotter_expressions(shifted)
        _merge_trotter_stage!(stages, key, expression)
    end
    for (sector, coefficient) in grouped_trotter_constants(shifted)
        key = TrotterFrequencyGroupKey(sector, Int[])
        _merge_trotter_stage!(
            stages, key, _constant_trotter_expression(coefficient, sector)
        )
    end

    blocked = Dict{TrotterFrequencyGroupKey{S},TrotterFrequencyBlockedContribution{D,S}}()

    while !isempty(stages)
        key = _maximal_trotter_stage(stages)
        active = key.active_frequency_indices
        isempty(active) && break

        expression = pop!(stages, key)
        frequency_index = _next_trotter_frequency(stages, key)
        plan = CausalFrequencyReductionPlan((frequency_index,))
        result = reduce_causal_frequency_with_support(expression, plan)
        kind = causal_frequency_support_kind(result)

        if kind === CausalFrequencyIntegrated
            next_active = Int[index for index in active if index != frequency_index]
            next_key = TrotterFrequencyGroupKey(
                canonical_frequency_sector(key), next_active
            )
            _merge_trotter_stage!(
                stages, next_key, causal_frequency_support_expression(result)
            )
            continue
        end

        blocked[key] = TrotterFrequencyBlockedContribution{D,S}(
            key,
            causal_frequency_support_expression(result),
            kind,
            causal_frequency_support_blocked_frequency(result),
            causal_frequency_singular_supports(result),
        )
    end

    constants = Dict{CanonicalFrequencySector{S},D}()
    boundary = Dict{CanonicalFrequencySector{S},CausalFrequencyExpression{D,S}}()
    for (key, expression) in stages
        isempty(key.active_frequency_indices) ||
            throw(ErrorException("active causal stage escaped canonical Trotter reduction"))
        sector = canonical_frequency_sector(key)
        if _is_trotter_constant_expression(expression)
            coefficient = _trotter_constant_coefficient(expression)
            iszero(coefficient) || (constants[sector] = coefficient)
        else
            boundary[sector] = expression
        end
    end

    return CanonicalTrotterFrequencyReduction{D,S,E1,E2}(
        constants, boundary, blocked, copy(unresolved_trotter_states(shifted))
    )
end
