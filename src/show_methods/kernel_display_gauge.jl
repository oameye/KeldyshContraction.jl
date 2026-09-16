# Display-only reconstruction of a coherent dummy-loop gauge for quotient kernels.
#
# `quotient_loop_momenta` intentionally canonicalizes each occupation-monomial ×
# kinematic-monomial atom independently. That is the correct canonical IR, but different atoms
# can therefore use different representatives of the same signed-permutation orbit. Before
# looking for a physics-facing factorization of the final kernel, first keep any shorter exact
# factorization already present in the stored quotient gauge. Only if that direct representation
# cannot be compacted do we lift quotient atoms into one common support gauge and average over
# the residual dummy-loop stabilizer. This never mutates the compiler IR and uses only the same
# universally safe signed-permutation subgroup as the production quotient.

function _display_loop_transforms(basis::MomentumBasis, external::MomentumVariable)
    nloops = length(basis) - 1
    nloops <= 4 || return LoopMomentumTransform[]

    transforms = LoopMomentumTransform[]
    slots = collect(1:nloops)
    external_shifts = zeros(Int, nloops)
    for permutation in Combinatorics.permutations(slots)
        permutation_vector = Int[value for value in permutation]
        for mask in 0:((1 << nloops) - 1)
            signs = Vector{Int}(undef, nloops)
            for slot in 1:nloops
                signs[slot] = iszero(mask & (1 << (slot - 1))) ? 1 : -1
            end
            push!(
                transforms,
                loop_permutation_transform(
                    basis, external, permutation_vector, signs, external_shifts
                ),
            )
        end
    end
    return transforms
end

function _display_parameter_isless(a::ParameterMonomial, b::ParameterMonomial)
    n = min(length(a.powers), length(b.powers))
    @inbounds for i in 1:n
        ap = a.powers[i]
        bp = b.powers[i]
        ap.name === bp.name || return isless(ap.name, bp.name)
        ap.exponent == bp.exponent || return ap.exponent < bp.exponent
    end
    return length(a.powers) < length(b.powers)
end

function _display_basis_isless(a::MomentumBasis, b::MomentumBasis)
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        a[i] == b[i] || return isless(a[i], b[i])
    end
    return length(a) < length(b)
end

function _display_number_isless(a::Number, b::Number)
    isequal(a, b) && return false
    ar = real(a)
    br = real(b)
    isequal(ar, br) || return isless(ar, br)
    return isless(imag(a), imag(b))
end

function _display_support_isless(a::FrequencySupport{S}, b::FrequencySupport{S}) where {S}
    nshells = min(length(a.shells), length(b.shells))
    @inbounds for i in 1:nshells
        ashell = a.shells[i]
        bshell = b.shells[i]
        isequal(ashell, bshell) || return isless(ashell, bshell)
    end
    length(a.shells) == length(b.shells) || return length(a.shells) < length(b.shells)

    npv = min(length(a.principal_values), length(b.principal_values))
    @inbounds for i in 1:npv
        apv = a.principal_values[i]
        bpv = b.principal_values[i]
        isequal(apv, bpv) || return isless(apv, bpv)
    end
    return length(a.principal_values) < length(b.principal_values)
end

function _display_kinematic_isless(a::MomentumPolynomial, b::MomentumPolynomial)
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        am, ac = a.terms[i]
        bm, bc = b.terms[i]
        isequal(am, bm) || return isless(am, bm)
        isequal(ac, bc) || return _display_number_isless(ac, bc)
    end
    return length(a) < length(b)
end

function _display_occupation_isless(
    a::OccupationPolynomial{C1,S}, b::OccupationPolynomial{C2,S}
) where {C1<:Number,C2<:Number,S<:Statistics}
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        am, ac = a.terms[i]
        bm, bc = b.terms[i]
        isequal(am, bm) || return isless(am, bm)
        isequal(ac, bc) || return _display_number_isless(ac, bc)
    end
    return length(a) < length(b)
end

function _kernel_display_pair_isless(
    a::Pair{CollisionKernelSector{S},OccupationPolynomial{C1,S}},
    b::Pair{CollisionKernelSector{S},OccupationPolynomial{C2,S}},
) where {C1<:Number,C2<:Number,S<:Statistics}
    a_sector, a_polynomial = a
    b_sector, b_polynomial = b

    a_parameter = parameters(a_sector)
    b_parameter = parameters(b_sector)
    isequal(a_parameter, b_parameter) ||
        return _display_parameter_isless(a_parameter, b_parameter)

    a_basis = momentum_basis(a_sector)
    b_basis = momentum_basis(b_sector)
    isequal(a_basis, b_basis) || return _display_basis_isless(a_basis, b_basis)

    a_external = external_wigner_momentum(a_sector)
    b_external = external_wigner_momentum(b_sector)
    a_external == b_external || return isless(a_external, b_external)

    a_support = frequency_support(a_sector)
    b_support = frequency_support(b_sector)
    isequal(a_support, b_support) || return _display_support_isless(a_support, b_support)

    a_kinematic = kinematic_factor(a_sector)
    b_kinematic = kinematic_factor(b_sector)
    isequal(a_kinematic, b_kinematic) ||
        return _display_kinematic_isless(a_kinematic, b_kinematic)

    return _display_occupation_isless(a_polynomial, b_polynomial)
end

struct _DisplaySupportOrbit{S<:Statistics}
    canonical::FrequencySupport{S}
    images::Dict{
        FrequencySupport{S},Vector{Tuple{LoopMomentumTransform,MomentumCoefficient}}
    }
end

function _display_support_orbit(
    support::FrequencySupport{S}, transforms::Vector{LoopMomentumTransform}
) where {S<:Statistics}
    images = Dict{
        FrequencySupport{S},Vector{Tuple{LoopMomentumTransform,MomentumCoefficient}}
    }()
    if isempty(transforms)
        return _DisplaySupportOrbit{S}(support, images)
    end

    canonical = support
    for transform in transforms
        transformed, factor = _transform_frequency_support(support, transform)
        matches = get!(images, transformed) do
            return Tuple{LoopMomentumTransform,MomentumCoefficient}[]
        end
        push!(matches, (transform, factor))
        _display_support_isless(transformed, canonical) && (canonical = transformed)
    end
    return _DisplaySupportOrbit{S}(canonical, images)
end

function _display_support_orbit_cache(
    group::Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}},
    transforms::Vector{LoopMomentumTransform},
) where {C<:Number,S<:Statistics}
    cache = Dict{FrequencySupport{S},_DisplaySupportOrbit{S}}()
    for (sector, _) in group
        support = frequency_support(sector)
        haskey(cache, support) && continue
        cache[support] = _display_support_orbit(support, transforms)
    end
    return cache
end

function _display_support_orbit_groups(
    group::Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}},
    cache::Dict{FrequencySupport{S},_DisplaySupportOrbit{S}},
) where {C<:Number,S<:Statistics}
    groups = Dict{
        FrequencySupport{S},Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}}
    }()
    for pair in group
        support = frequency_support(first(pair))
        canonical = cache[support].canonical
        orbit_group = get!(groups, canonical) do
            return Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}[]
        end
        push!(orbit_group, pair)
    end

    supports = collect(keys(groups))
    sort!(supports; lt=_display_support_isless)
    orbit_groups = Vector{Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}}}()
    sizehint!(orbit_groups, length(supports))
    for support in supports
        orbit_group = groups[support]
        sort!(orbit_group; lt=_kernel_display_pair_isless)
        push!(orbit_groups, orbit_group)
    end
    return orbit_groups
end

function _lift_kernel_group_to_support(
    group::Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}},
    target_support::FrequencySupport{S},
    cache::Dict{FrequencySupport{S},_DisplaySupportOrbit{S}},
) where {C<:Number,S<:Statistics}
    rows = Dict{CollisionKernelSector{S},OccupationPolynomial{ComplexRationals,S}}()
    empty_rows = Pair{CollisionKernelSector{S},OccupationPolynomial{ComplexRationals,S}}[]

    for (sector, polynomial) in group
        orbit = cache[frequency_support(sector)]
        matches = get(orbit.images, target_support, nothing)
        matches === nothing && return false, empty_rows
        orbit_weight = inv(convert(ComplexRationals, length(matches)))

        exact_terms = Pair{OccupationMonomial{S},ComplexRationals}[]
        sizehint!(exact_terms, length(polynomial))
        for (occupation_monomial, occupation_coefficient) in polynomial
            valid_coefficient, exact_coefficient = _complex_rational(occupation_coefficient)
            valid_coefficient || return false, empty_rows
            push!(exact_terms, occupation_monomial => exact_coefficient)
        end

        for (transform, support_factor) in matches
            transformed_kinematic = transform_loop_momenta(
                kinematic_factor(sector), transform
            )
            transform_factor =
                convert(ComplexRationals, support_factor) * orbit_weight

            for (occupation_monomial, exact_coefficient) in exact_terms
                transformed_occupation = transform_loop_momenta(
                    occupation_monomial, transform
                )
                for (kinematic_monomial, kinematic_coefficient) in transformed_kinematic
                    unit_kinematic = MomentumPolynomial(
                        kinematic_monomial, one(ComplexRationals)
                    )
                    row_sector = CollisionKernelSector{S}(
                        parameters(sector),
                        momentum_basis(sector),
                        external_wigner_momentum(sector),
                        unit_kinematic,
                        target_support,
                    )
                    coefficient =
                        exact_coefficient * transform_factor * kinematic_coefficient
                    contribution = OccupationPolynomial{ComplexRationals,S}(
                        Pair{OccupationMonomial{S},ComplexRationals}[
                            transformed_occupation => coefficient
                        ],
                    )
                    if haskey(rows, row_sector)
                        combined = rows[row_sector] + contribution
                        if iszero(combined)
                            delete!(rows, row_sector)
                        else
                            rows[row_sector] = combined
                        end
                    elseif !iszero(contribution)
                        rows[row_sector] = contribution
                    end
                end
            end
        end
    end

    lifted = collect(rows)
    sort!(lifted; lt=_kernel_display_pair_isless)
    return true, lifted
end

function _exact_kernel_rank_one_kinematic(
    group::Vector{
        Pair{CollisionKernelSector{S},OccupationPolynomial{ComplexRationals,S}}
    },
) where {S<:Statistics}
    reference = last(first(group))
    kinematic = MomentumPolynomial{ComplexRationals}()
    for (sector, polynomial) in group
        valid_scale, scale = _polynomial_scale(polynomial, reference)
        valid_scale || return false, kinematic, reference
        kinematic += scale * kinematic_factor(sector)
    end
    return true, kinematic, reference
end

function _aligned_rank_one_kinematic(
    group::Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}},
    cache::Dict{FrequencySupport{S},_DisplaySupportOrbit{S}},
)::Tuple{
    Bool,
    MomentumPolynomial{ComplexRationals},
    OccupationPolynomial{ComplexRationals,S},
    FrequencySupport{S},
} where {C<:Number,S<:Statistics}
    first_sector = first(first(group))
    target_support = frequency_support(first_sector)
    empty_kinematic = MomentumPolynomial{ComplexRationals}()
    empty_distribution = OccupationPolynomial{ComplexRationals,S}()

    lifted_ok, lifted = _lift_kernel_group_to_support(group, target_support, cache)
    lifted_ok || return false, empty_kinematic, empty_distribution, target_support
    isempty(lifted) && return false, empty_kinematic, empty_distribution, target_support

    rank_one, kinematic, distribution = _exact_kernel_rank_one_kinematic(lifted)
    return rank_one, kinematic, distribution, target_support
end

function _rank_one_display_candidate(
    sector::CollisionKernelSector,
    kinematic::MomentumPolynomial{ComplexRationals},
    distribution::OccupationPolynomial,
    support::FrequencySupport,
    latex::Bool,
)
    basis = momentum_basis(sector)
    external = external_wigner_momentum(sector)
    coefficient, kinematic_text = _compact_kinematic_string(
        kinematic, distribution, support, basis, external, latex
    )
    distribution_text = _compact_distribution_string(distribution, basis, external, latex)
    factors = String[
        _parameter_string(parameters(sector), latex),
        _measure_string(basis, external, latex),
        kinematic_text,
        distribution_text == "1" ? "" : distribution_text,
        _frequency_support_string(support, basis, external, latex),
    ]
    return _term_string(coefficient, factors, latex)
end

function _compact_group_string(
    group::Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}},
    latex::Bool,
    cache::Dict{FrequencySupport{S},_DisplaySupportOrbit{S}},
) where {C<:Number,S<:Statistics}
    isempty(group) && return "0"
    sector = first(first(group))
    support = frequency_support(sector)
    expanded = _physical_collision_string(group, latex)

    # Preserve a compact exact expression already visible in the quotient gauge. This is
    # essential for support-free one-loop kernels, where averaging over the full support
    # stabilizer would erase the physically useful routed-momentum representative.
    literal_support = all(pair -> frequency_support(first(pair)) == support, group)
    if literal_support
        rank_one, kinematic, distribution = _rank_one_kinematic(group)
        if rank_one
            direct = _rank_one_display_candidate(
                sector, kinematic, distribution, support, latex
            )
            _display_cost(direct) < _display_cost(expanded) && return direct
        end
    end

    # The direct quotient gauge was not compact enough. Reconstruct a coherent gauge across the
    # signed-permutation orbit and retry the same exact rank-one/factorization machinery.
    rank_one, kinematic, distribution, aligned_support = _aligned_rank_one_kinematic(
        group, cache
    )
    rank_one || return expanded
    compact = _rank_one_display_candidate(
        sector, kinematic, distribution, aligned_support, latex
    )
    return _display_cost(compact) < _display_cost(expanded) ? compact : expanded
end

# Kernel sectors must be classified by support orbit before compaction. The generic display
# path groups by literal support, which is appropriate before the loop quotient but would split
# physically identical quotient representatives before their common gauge can be reconstructed.
function _compact_physical_collision_string_impl(
    terms::Dict{CollisionKernelSector{S},OccupationPolynomial{C,S}}, ::Type{S}, latex::Bool
) where {C<:Number,S<:Statistics}
    BaseKey = Tuple{ParameterMonomial,MomentumBasis,MomentumVariable}
    base_groups = Dict{
        BaseKey,Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}}
    }()

    for (sector, polynomial) in terms
        key = (parameters(sector), momentum_basis(sector), external_wigner_momentum(sector))
        group = get!(base_groups, key) do
            return Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}[]
        end
        push!(group, sector => polynomial)
    end

    rendered = String[]
    for base_group in values(base_groups)
        sector = first(first(base_group))
        transforms = _display_loop_transforms(
            momentum_basis(sector), external_wigner_momentum(sector)
        )
        cache = _display_support_orbit_cache(base_group, transforms)
        for orbit_group in _display_support_orbit_groups(base_group, cache)
            push!(rendered, _compact_group_string(orbit_group, latex, cache))
        end
    end
    sort!(rendered)
    return _sum_string(rendered, latex)
end
