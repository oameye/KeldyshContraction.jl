# Display-only reconstruction of a coherent dummy-loop gauge for quotient kernels.
#
# `quotient_loop_momenta` intentionally canonicalizes each complete integrand atom
# independently.  That is the correct canonical IR, but different atoms can therefore use
# different representatives of the same signed-permutation orbit.  Before looking for a
# physics-facing factorization of the final kernel, align those representatives into one
# common loop gauge.  This never mutates the compiler IR and uses only the same universally
# safe signed-permutation subgroup as the production quotient.

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

function _display_transform_kernel_pair(
    sector::CollisionKernelSector{S},
    polynomial::OccupationPolynomial{C,S},
    transform::LoopMomentumTransform,
) where {C<:Number,S<:Statistics}
    support, support_factor = _transform_frequency_support(
        frequency_support(sector), transform
    )
    kinematic = transform_loop_momenta(kinematic_factor(sector), transform)
    distribution =
        convert(C, support_factor) * transform_loop_momenta(polynomial, transform)
    return kinematic, distribution, support
end

function _aligned_rank_one_kinematic(
    group::Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}}
) where {C<:Number,S<:Statistics}
    first_sector, first_polynomial = first(group)
    basis = momentum_basis(first_sector)
    external = external_wigner_momentum(first_sector)
    transforms = _display_loop_transforms(basis, external)

    empty_kinematic = MomentumPolynomial{ComplexRationals}()
    isempty(transforms) &&
        return false, empty_kinematic, first_polynomial, frequency_support(first_sector)

    for reference_transform in transforms
        reference_kinematic, reference_distribution, reference_support = _display_transform_kernel_pair(
            first_sector, first_polynomial, reference_transform
        )
        kinematic = reference_kinematic
        aligned = true

        for index in 2:length(group)
            sector, polynomial = group[index]
            matched = false
            for transform in transforms
                candidate_kinematic, candidate_distribution, candidate_support = _display_transform_kernel_pair(
                    sector, polynomial, transform
                )
                candidate_support == reference_support || continue

                valid_scale, scale = _polynomial_scale(
                    candidate_distribution, reference_distribution
                )
                valid_scale || continue
                valid_coefficient, coefficient = _complex_rational(scale)
                valid_coefficient || continue

                kinematic += coefficient * candidate_kinematic
                matched = true
                break
            end
            if !matched
                aligned = false
                break
            end
        end

        aligned && return true, kinematic, reference_distribution, reference_support
    end

    return false, empty_kinematic, first_polynomial, frequency_support(first_sector)
end

function _compact_group_string(
    group::Vector{Pair{CollisionKernelSector{S},OccupationPolynomial{C,S}}}, latex::Bool
) where {C<:Number,S<:Statistics}
    isempty(group) && return "0"
    sector = first(first(group))
    basis = momentum_basis(sector)
    external = external_wigner_momentum(sector)

    rank_one, kinematic, distribution, support = _aligned_rank_one_kinematic(group)
    rank_one || return _physical_collision_string(group, latex)

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
    compact = _term_string(coefficient, factors, latex)
    expanded = _physical_collision_string(group, latex)
    return _display_cost(compact) < _display_cost(expanded) ? compact : expanded
end
