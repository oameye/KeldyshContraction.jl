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
                    basis, external, permutation_vector, signs, zeros(Int, nloops)
                ),
            )
        end
    end
    return transforms
end

function _kernel_display_pair_key(
    pair::Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}
) where {C<:Number,S<:Statistics}
    singleton = Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}[pair]
    return _physical_collision_string(singleton, false)
end

function _display_support_equivalent(
    source::FrequencySupport{S},
    target::FrequencySupport{S},
    transforms::Vector{LoopMomentumTransform},
) where {S<:Statistics}
    isempty(transforms) && return source == target
    for transform in transforms
        transformed, _ = _transform_frequency_support(source, transform)
        transformed == target && return true
    end
    return false
end

function _display_support_orbit_groups(
    group::Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}},
    transforms::Vector{LoopMomentumTransform},
) where {C<:Number,S<:Statistics}
    ordered = sort!(copy(group); by=_kernel_display_pair_key)
    orbit_groups = Vector{Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}}}()

    for pair in ordered
        support = frequency_support(first(pair))
        placed = false
        for orbit_group in orbit_groups
            reference_support = frequency_support(first(first(orbit_group)))
            _display_support_equivalent(support, reference_support, transforms) || continue
            push!(orbit_group, pair)
            placed = true
            break
        end
        if !placed
            push!(
                orbit_groups, Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}[pair]
            )
        end
    end
    return orbit_groups
end

function _matching_support_transforms(
    support::FrequencySupport{S},
    target::FrequencySupport{S},
    transforms::Vector{LoopMomentumTransform},
) where {S<:Statistics}
    matches = Tuple{LoopMomentumTransform,MomentumCoefficient}[]
    for transform in transforms
        transformed, factor = _transform_frequency_support(support, transform)
        transformed == target || continue
        push!(matches, (transform, factor))
    end
    return matches
end

function _lift_kernel_group_to_support(
    group::Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}},
    target_support::FrequencySupport{S},
    transforms::Vector{LoopMomentumTransform},
) where {C<:Number,S<:Statistics}
    rows = Dict{CollisionKernelSector{S},OccupationPolynomial{ComplexRationals,S}}()
    empty_rows = Pair{
        CollisionKernelSector{S},OccupationPolynomial{ComplexRationals,S}
    }[]

    for (sector, polynomial) in group
        matches = _matching_support_transforms(
            frequency_support(sector), target_support, transforms
        )
        isempty(matches) && return false, empty_rows
        orbit_weight = inv(convert(ComplexRationals, length(matches)))

        for (occupation_monomial, occupation_coefficient) in polynomial
            valid_coefficient, exact_coefficient = _complex_rational(occupation_coefficient)
            valid_coefficient || return false, empty_rows

            for (transform, support_factor) in matches
                transformed_occupation = transform_loop_momenta(
                    occupation_monomial, transform
                )
                transformed_kinematic = transform_loop_momenta(
                    kinematic_factor(sector), transform
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
                        exact_coefficient *
                        convert(ComplexRationals, support_factor) *
                        kinematic_coefficient *
                        orbit_weight
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
    sort!(lifted; by=_kernel_display_pair_key)
    return true, lifted
end

function _aligned_rank_one_kinematic(
    group::Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}}
) where {C<:Number,S<:Statistics}
    first_sector = first(first(group))
    basis = momentum_basis(first_sector)
    external = external_wigner_momentum(first_sector)
    transforms = _display_loop_transforms(basis, external)
    empty_kinematic = MomentumPolynomial{ComplexRationals}()

    isempty(transforms) &&
        return false, empty_kinematic, last(first(group)), frequency_support(first_sector)

    ordered = sort!(copy(group); by=_kernel_display_pair_key)
    target_support = frequency_support(first(first(ordered)))
    lifted_ok, lifted = _lift_kernel_group_to_support(ordered, target_support, transforms)
    lifted_ok || return false, empty_kinematic, last(first(ordered)), target_support
    isempty(lifted) && return false, empty_kinematic, last(first(ordered)), target_support

    rank_one, kinematic, distribution = _rank_one_kinematic(lifted)
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
    group::Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}}, latex::Bool
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
    rank_one, kinematic, distribution, aligned_support = _aligned_rank_one_kinematic(group)
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
        ordered = sort!(copy(base_group); by=_kernel_display_pair_key)
        sector = first(first(ordered))
        transforms = _display_loop_transforms(
            momentum_basis(sector), external_wigner_momentum(sector)
        )
        for orbit_group in _display_support_orbit_groups(ordered, transforms)
            push!(rendered, _compact_group_string(orbit_group, latex))
        end
    end
    sort!(rendered)
    return _sum_string(rendered, latex)
end
