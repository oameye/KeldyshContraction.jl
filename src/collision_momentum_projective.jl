# Projective normalization of routed derivative factors at the physical loop quotient.
#
# A MomentumComponent stores an exact linear form such as q or -q. Scalar multiples
# of the same linear form are algebraically the same polynomial generator; the scalar
# belongs in the MomentumPolynomial coefficient. Canonicalizing that scalar before a
# dummy-loop gauge is chosen prevents coordinate orientation from erasing derivative
# signs.

function _projective_kinematic_monomial(monomial::MomentumMonomial)
    components = MomentumComponent[]
    sizehint!(components, length(monomial))
    factor = one(MomentumCoefficient)

    for component in monomial
        pivot = zero(MomentumCoefficient)
        for coefficient in component.momentum.coefficients
            if !iszero(coefficient)
                pivot = coefficient
                break
            end
        end
        iszero(pivot) && return MomentumMonomial(), zero(MomentumCoefficient), false

        normalized_momentum = inv(pivot) * component.momentum
        push!(components, MomentumComponent(normalized_momentum, component.axis))
        factor *= pivot
    end

    return MomentumMonomial(components), factor, true
end

# Derivative lowering uses this exact coefficient domain throughout the kinetic path.
# Put each routed linear factor in projective normal form at construction time so the
# later loop quotient sees q_x and -q_x as one generator with coefficients +1 and -1,
# rather than as two coordinate spellings whose signs can be absorbed by different
# dummy-loop gauges.
function MomentumPolynomial(
    monomial::MomentumMonomial, coefficient::ComplexRationals
)::MomentumPolynomial{ComplexRationals}
    normalized, factor, nonzero = _projective_kinematic_monomial(monomial)
    (nonzero && !iszero(coefficient)) || return MomentumPolynomial{ComplexRationals}()
    normalized_coefficient = coefficient * convert(ComplexRationals, factor)
    return MomentumPolynomial{ComplexRationals}([normalized => normalized_coefficient])
end

function transform_loop_momenta(
    polynomial::MomentumPolynomial{ComplexRationals}, transform::LoopMomentumTransform
)
    terms = Pair{MomentumMonomial,ComplexRationals}[]
    sizehint!(terms, length(polynomial))

    for (monomial, coefficient) in polynomial
        transformed = transform_loop_momenta(monomial, transform)
        normalized, factor, nonzero = _projective_kinematic_monomial(transformed)
        nonzero || continue
        push!(terms, normalized => coefficient * convert(ComplexRationals, factor))
    end

    return MomentumPolynomial{ComplexRationals}(terms)
end

# Return the primitive integer representative of one routed linear form without
# choosing between L and -L. Magnitudes are graph semantics; the global orientation
# is chosen only by the complete loop-canonicalization witness.
function _primitive_projective_momentum(momentum::LinearMomentum)
    common_denominator = 1
    has_nonzero = false
    for coefficient in momentum.coefficients
        iszero(coefficient) && continue
        has_nonzero = true
        common_denominator = lcm(common_denominator, denominator(coefficient))
    end
    has_nonzero || return zeros(MomentumCoefficient, length(momentum.coefficients))

    integers = Vector{Int}(undef, length(momentum.coefficients))
    divisor = 0
    for index in eachindex(momentum.coefficients)
        coefficient = momentum.coefficients[index]
        value = numerator(coefficient) * (common_denominator ÷ denominator(coefficient))
        integers[index] = value
        iszero(value) || (divisor = gcd(divisor, abs(value)))
    end

    primitive = Vector{MomentumCoefficient}(undef, length(integers))
    for index in eachindex(integers)
        primitive[index] = (integers[index] ÷ divisor) // 1
    end
    return primitive
end

# Encode a derivative linear form as the unoriented projective line {L,-L}. The two
# local orientation vertices have identical colors; coefficient magnitudes connect
# them to the fixed external orientation or to a dummy loop's two global orientations.
# A signed loop relabeling therefore acts on this graph covariantly without turning
# the physical scalar coefficient of the kinetic term into graph identity.
function _add_projective_kinematic_component!(
    builder::_LoopCanonicalGraphBuilder,
    term_vertex::Int,
    component::MomentumComponent,
    external_index::Int,
    loop_indices::Vector{Int},
    positive_vertices::Vector{Int},
    negative_vertices::Vector{Int},
    external_positive::Int,
    external_negative::Int,
)
    component_pair = _add_loop_vertex!(builder, _loop_axis_color(21, component.axis))
    component_positive = _add_loop_vertex!(builder, _loop_graph_color(22))
    component_negative = _add_loop_vertex!(builder, _loop_graph_color(22))
    _add_loop_edge!(builder, term_vertex, component_pair)
    _add_loop_edge!(builder, component_pair, component_positive)
    _add_loop_edge!(builder, component_pair, component_negative)

    primitive = _primitive_projective_momentum(component.momentum)
    for basis_index in eachindex(primitive)
        coefficient = primitive[basis_index]
        iszero(coefficient) && continue
        magnitude = abs(coefficient)

        loop_slot = 0
        if basis_index != external_index
            for (slot, loop_index) in enumerate(loop_indices)
                if basis_index == loop_index
                    loop_slot = slot
                    break
                end
            end
            iszero(loop_slot) && error("projective derivative uses an unknown momentum basis slot")
        end

        for local_sign in (1, -1)
            incidence = _add_loop_vertex!(
                builder,
                _loop_axis_color(23, :projective_component_coefficient, magnitude),
            )
            component_orientation = local_sign == 1 ? component_positive : component_negative
            coefficient_positive = (local_sign == 1) == (coefficient > 0)
            target = if basis_index == external_index
                coefficient_positive ? external_positive : external_negative
            else
                coefficient_positive ? positive_vertices[loop_slot] : negative_vertices[loop_slot]
            end
            _add_loop_edge!(builder, component_orientation, incidence)
            _add_loop_edge!(builder, incidence, target)
        end
    end
    return builder
end

function _add_projective_loop_semantics!(
    builder::_LoopCanonicalGraphBuilder,
    root::Int,
    sector::ReducedCollisionSector{S},
    monomial::OccupationMonomial{S},
    external_index::Int,
    loop_indices::Vector{Int},
    positive_vertices::Vector{Int},
    negative_vertices::Vector{Int},
    loop_incidences::Vector{Vector{Tuple{Int,MomentumCoefficient}}},
) where {S<:Statistics}
    external_positive = _add_loop_vertex!(builder, _loop_graph_color(5))
    external_negative = _add_loop_vertex!(builder, _loop_graph_color(6))
    _add_loop_edge!(builder, root, external_positive)
    _add_loop_edge!(builder, root, external_negative)

    for atom in monomial
        momentum = atom.momentum
        vertex = _add_loop_vertex!(
            builder, _loop_family_color(10, atom.family, momentum[external_index])
        )
        _add_loop_edge!(builder, root, vertex)
        _add_loop_momentum_incidence!(
            builder,
            vertex,
            momentum,
            loop_indices,
            positive_vertices,
            negative_vertices,
            loop_incidences,
        )
    end

    for (kinematic_monomial, _) in kinematic_factor(sector)
        # The physical scalar coefficient is a linear weight, not canonical graph data.
        term_vertex = _add_loop_vertex!(builder, _loop_graph_color(20))
        _add_loop_edge!(builder, root, term_vertex)
        for component in kinematic_monomial
            _add_projective_kinematic_component!(
                builder,
                term_vertex,
                component,
                external_index,
                loop_indices,
                positive_vertices,
                negative_vertices,
                external_positive,
                external_negative,
            )
        end
    end

    support = frequency_support(sector)
    for shell in support.shells
        support_vertex = _add_loop_vertex!(builder, _loop_graph_color(30))
        _add_loop_edge!(builder, root, support_vertex)
        for (atom, coefficient) in shell.energy
            momentum = atom.momentum
            energy_vertex = _add_loop_vertex!(
                builder,
                _loop_family_color(32, atom.family, coefficient, momentum[external_index]),
            )
            _add_loop_edge!(builder, support_vertex, energy_vertex)
            _add_loop_momentum_incidence!(
                builder,
                energy_vertex,
                momentum,
                loop_indices,
                positive_vertices,
                negative_vertices,
                loop_incidences,
            )
        end
    end
    for principal_value in support.principal_values
        support_vertex = _add_loop_vertex!(builder, _loop_graph_color(31))
        _add_loop_edge!(builder, root, support_vertex)
        for (atom, coefficient) in principal_value.energy
            momentum = atom.momentum
            energy_vertex = _add_loop_vertex!(
                builder,
                _loop_family_color(32, atom.family, coefficient, momentum[external_index]),
            )
            _add_loop_edge!(builder, support_vertex, energy_vertex)
            _add_loop_momentum_incidence!(
                builder,
                energy_vertex,
                momentum,
                loop_indices,
                positive_vertices,
                negative_vertices,
                loop_incidences,
            )
        end
    end
    return builder
end

function _projective_canonical_loop_transform(
    sector::ReducedCollisionSector{S}, monomial::OccupationMonomial{S}
) where {S<:Statistics}
    basis = momentum_basis(sector)
    external = external_wigner_momentum(sector)
    external_index = _external_basis_index(basis, external)
    loop_indices = _loop_basis_indices(basis, external)
    nloops = length(loop_indices)

    builder = _LoopCanonicalGraphBuilder()
    root = _add_loop_vertex!(builder, _loop_graph_color(1))
    pair_vertices = Vector{Int}(undef, nloops)
    positive_vertices = Vector{Int}(undef, nloops)
    negative_vertices = Vector{Int}(undef, nloops)
    loop_incidences = [Tuple{Int,MomentumCoefficient}[] for _ in 1:nloops]
    for slot in 1:nloops
        pair = _add_loop_vertex!(builder, _loop_graph_color(2))
        positive = _add_loop_vertex!(builder, _loop_graph_color(3))
        negative = _add_loop_vertex!(builder, _loop_graph_color(3))
        pair_vertices[slot] = pair
        positive_vertices[slot] = positive
        negative_vertices[slot] = negative
        _add_loop_edge!(builder, root, pair)
        _add_loop_edge!(builder, pair, positive)
        _add_loop_edge!(builder, pair, negative)
    end

    _add_projective_loop_semantics!(
        builder,
        root,
        sector,
        monomial,
        external_index,
        loop_indices,
        positive_vertices,
        negative_vertices,
        loop_incidences,
    )

    permutation = NautyGraphs.canonical_permutation(_loop_graph(builder))
    rank = Vector{Int}(undef, length(permutation))
    for (canonical_rank, original_vertex) in enumerate(permutation)
        rank[original_vertex] = canonical_rank
    end

    ordered_slots = sortperm(1:nloops; by=slot -> rank[pair_vertices[slot]])
    loop_permutation = zeros(Int, nloops)
    loop_signs = ones(Int, nloops)
    for (new_slot, old_slot) in enumerate(ordered_slots)
        loop_permutation[old_slot] = new_slot
        loop_signs[old_slot] =
            rank[positive_vertices[old_slot]] < rank[negative_vertices[old_slot]] ? 1 : -1
    end
    return loop_permutation_transform(
        basis, external, loop_permutation, loop_signs, zeros(Int, nloops)
    )
end

# These are the two physical statistics supported by the package. Keeping the projective
# semantic override here avoids changing the generic graph machinery: both current and future
# canonical-labeling backends consume the same KC-owned semantic graph.
function _canonical_loop_transform(
    sector::ReducedCollisionSector{Boson}, monomial::OccupationMonomial{Boson}
)
    return _projective_canonical_loop_transform(sector, monomial)
end
function _canonical_loop_transform(
    sector::ReducedCollisionSector{Fermion}, monomial::OccupationMonomial{Fermion}
)
    return _projective_canonical_loop_transform(sector, monomial)
end

# Signed dummy-loop reflections are part of the quotient group. If one such
# reflection leaves the occupation monomial and frequency support unchanged while
# multiplying the derivative kinematics by c != 1, the complete atom satisfies
# I = c I and therefore vanishes exactly. This does not assume n(q) = n(-q): atoms
# containing the reflected occupation momentum fail the invariance test and survive.
function _vanishes_under_loop_reflection(
    sector::CollisionKernelSector{S}, monomial::OccupationMonomial{S}
) where {S<:Statistics}
    basis = momentum_basis(sector)
    external = external_wigner_momentum(sector)
    nloops = length(basis) - 1
    iszero(nloops) && return false

    kinematic = kinematic_factor(sector)
    length(kinematic) == 1 || return false
    reference_monomial, reference_coefficient = only(kinematic)
    reference_coefficient == one(ComplexRationals) || return false

    permutation = collect(1:nloops)
    shifts = zeros(Int, nloops)
    signs = ones(Int, nloops)
    for slot in 1:nloops
        signs[slot] = -1
        transform = loop_permutation_transform(basis, external, permutation, signs, shifts)
        transformed_monomial = transform_loop_momenta(monomial, transform)
        if transformed_monomial == monomial
            transformed_support, support_factor = _transform_frequency_support(
                frequency_support(sector), transform
            )
            if transformed_support == frequency_support(sector)
                transformed_kinematic = transform_loop_momenta(kinematic, transform)
                if length(transformed_kinematic) == 1
                    reflected_monomial, kinematic_factor = only(transformed_kinematic)
                    if reflected_monomial == reference_monomial
                        factor =
                            kinematic_factor * convert(ComplexRationals, support_factor)
                        factor != one(ComplexRationals) && return true
                    end
                end
            end
        end
        signs[slot] = 1
    end
    return false
end

function _merge_kernel_polynomial!(
    terms::Dict{CollisionKernelSector{S},OccupationPolynomial{C,S}},
    sector::CollisionKernelSector{S},
    polynomial::OccupationPolynomial{C,S},
) where {C<:Number,S<:Statistics}
    iszero(polynomial) && return terms
    if haskey(terms, sector)
        combined = terms[sector] + polynomial
        if iszero(combined)
            delete!(terms, sector)
        else
            terms[sector] = combined
        end
    else
        terms[sector] = polynomial
    end
    return terms
end

# `quotient_loop_momenta` promotes its output coefficient domain with
# `ComplexRationals`, so this method is the canonical insertion path for quotient
# atoms. Any scalar generated by the chosen loop gauge is moved from the transformed
# kinematic polynomial into the occupation-polynomial coefficient before sector
# identity is compared.
function _push_kernel_polynomial!(
    terms::Dict{CollisionKernelSector{S},OccupationPolynomial{C,S}},
    sector::CollisionKernelSector{S},
    polynomial::OccupationPolynomial{C,S},
) where {C<:Complex,S<:Statistics}
    iszero(polynomial) && return terms

    for (kinematic_monomial, kinematic_coefficient) in kinematic_factor(sector)
        iszero(kinematic_coefficient) && continue
        unit_kinematic = MomentumPolynomial(kinematic_monomial, one(ComplexRationals))
        unit_sector = CollisionKernelSector{S}(
            parameters(sector),
            momentum_basis(sector),
            external_wigner_momentum(sector),
            unit_kinematic,
            frequency_support(sector),
        )
        scale = convert(C, kinematic_coefficient)
        contributions = Pair{OccupationMonomial{S},C}[]
        sizehint!(contributions, length(polynomial))
        for (monomial, coefficient) in polynomial
            _vanishes_under_loop_reflection(unit_sector, monomial) && continue
            push!(contributions, monomial => scale * coefficient)
        end
        isempty(contributions) && continue
        _merge_kernel_polynomial!(
            terms, unit_sector, OccupationPolynomial{C,S}(contributions)
        )
    end

    return terms
end
