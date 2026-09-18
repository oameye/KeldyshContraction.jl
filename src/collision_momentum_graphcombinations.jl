# Research GraphCombinations backend for the physical dummy-loop quotient.
#
# The Nauty projective transform in `collision_momentum_projective.jl` remains an
# independent private oracle. Public Boson/Fermion quotienting is specialized below
# to use the exact GraphCombinations directed-canonicalization workspace.

function _prepare_loop_gc_labels!(
    labels::Vector{Int},
    color_order::Vector{_LoopGraphColor},
    builder::_LoopCanonicalGraphBuilder,
)::Int
    num_vertices = length(builder.colors)
    length(labels) >= num_vertices ||
        throw(DimensionMismatch("GC label buffer capacity is too small"))
    length(color_order) >= num_vertices ||
        throw(DimensionMismatch("GC color buffer capacity is too small"))

    copyto!(color_order, 1, builder.colors, 1, num_vertices)
    @inbounds for index in 2:num_vertices
        color = color_order[index]
        position = index - 1
        while position >= 1 && isless(color, color_order[position])
            color_order[position + 1] = color_order[position]
            position -= 1
        end
        color_order[position + 1] = color
    end

    @inbounds for vertex in eachindex(builder.colors)
        color = builder.colors[vertex]
        lower = 1
        upper = num_vertices + 1
        while lower < upper
            middle = (lower + upper) >>> 1
            if middle <= num_vertices && isless(color_order[middle], color)
                lower = middle + 1
            else
                upper = middle
            end
        end
        labels[vertex] = lower
    end
    return num_vertices
end

@static if isdefined(GC, :DirectedSimpleCanonicalizationWorkspace)
    struct _LoopGCCanonicalizationStorage
        workspace::GC.DirectedSimpleCanonicalizationWorkspace
        buffer::GC.DirectedCanonicalizationBuffer
        graph::GC.DirectedGCGraphBuffer
        labels::Vector{Int}
        color_order::Vector{_LoopGraphColor}
    end

    function _LoopGCCanonicalizationStorage(capacity::Int)
        return _LoopGCCanonicalizationStorage(
            GC.DirectedSimpleCanonicalizationWorkspace(
                capacity; frontier_capacity=max(4096, 64 * capacity)
            ),
            GC.DirectedCanonicalizationBuffer(capacity),
            GC.DirectedGCGraphBuffer(capacity),
            Vector{Int}(undef, capacity),
            Vector{_LoopGraphColor}(undef, capacity),
        )
    end

    function _canonicalize_loop_builder!(
        storage::_LoopGCCanonicalizationStorage, builder::_LoopCanonicalGraphBuilder
    )
        num_vertices = _prepare_loop_gc_labels!(
            storage.labels, storage.color_order, builder
        )
        GC.load_directed_graph!(storage.graph, builder.edges, num_vertices)
        GC.canonicalize_directed_simple!(
            storage.buffer, storage.workspace, storage.graph, storage.labels
        )
        return storage.buffer
    end
elseif isdefined(GC, :DirectedGCGraphBuffer)
    struct _LoopGCCanonicalizationStorage
        workspace::GC.DirectedCanonicalizationWorkspace
        buffer::GC.DirectedCanonicalizationBuffer
        graph::GC.DirectedGCGraphBuffer
        labels::Vector{Int}
        color_order::Vector{_LoopGraphColor}
    end

    function _LoopGCCanonicalizationStorage(capacity::Int)
        return _LoopGCCanonicalizationStorage(
            GC.DirectedCanonicalizationWorkspace(capacity),
            GC.DirectedCanonicalizationBuffer(capacity),
            GC.DirectedGCGraphBuffer(capacity),
            Vector{Int}(undef, capacity),
            Vector{_LoopGraphColor}(undef, capacity),
        )
    end

    function _canonicalize_loop_builder!(
        storage::_LoopGCCanonicalizationStorage, builder::_LoopCanonicalGraphBuilder
    )
        num_vertices = _prepare_loop_gc_labels!(
            storage.labels, storage.color_order, builder
        )
        GC.load_directed_graph!(storage.graph, builder.edges, num_vertices)
        GC.canonicalize_directed!(
            storage.buffer, storage.workspace, storage.graph, storage.labels
        )
        return storage.buffer
    end
else
    # Temporary research fallback for registered GraphCombinations v0.4.0. This branch
    # disappears once the capacity API from GraphCombinations#163 is released.
    struct _LoopGCCanonicalizationBuffers
        workspace::GC.DirectedCanonicalizationWorkspace
        buffer::GC.DirectedCanonicalizationBuffer
        graph::GC.DirectedGCGraph
        labels::Vector{Int}
        color_order::Vector{_LoopGraphColor}
    end

    function _LoopGCCanonicalizationBuffers(num_vertices::Int)
        return _LoopGCCanonicalizationBuffers(
            GC.DirectedCanonicalizationWorkspace(num_vertices),
            GC.DirectedCanonicalizationBuffer(num_vertices),
            GC.DirectedGCGraph(Pair{Int,Int}[], num_vertices),
            Vector{Int}(undef, num_vertices),
            Vector{_LoopGraphColor}(undef, num_vertices),
        )
    end

    struct _LoopGCCanonicalizationStorage
        cache::Dict{Int,_LoopGCCanonicalizationBuffers}
    end

    _LoopGCCanonicalizationStorage(::Int) =
        _LoopGCCanonicalizationStorage(Dict{Int,_LoopGCCanonicalizationBuffers}())

    function _canonicalize_loop_builder!(
        storage::_LoopGCCanonicalizationStorage, builder::_LoopCanonicalGraphBuilder
    )
        num_vertices = length(builder.colors)
        buffers = get!(storage.cache, num_vertices) do
            _LoopGCCanonicalizationBuffers(num_vertices)
        end
        _prepare_loop_gc_labels!(buffers.labels, buffers.color_order, builder)

        fill!(buffers.graph.multiplicities, 0)
        @inbounds for (source, target) in builder.edges
            slot = (source - 1) * num_vertices + target
            buffers.graph.multiplicities[slot] += 1
        end
        GC.canonicalize_directed!(
            buffers.buffer, buffers.workspace, buffers.graph, buffers.labels
        )
        return buffers.buffer
    end
end

mutable struct _LoopGCQuotientWorkspace
    canonicalization::_LoopGCCanonicalizationStorage
    builder::_LoopCanonicalGraphBuilder
    loop_indices::Vector{Int}
    pair_vertices::Vector{Int}
    positive_vertices::Vector{Int}
    negative_vertices::Vector{Int}
    loop_incidences::Vector{Vector{Tuple{Int,MomentumCoefficient}}}
    ordered_slots::Vector{Int}
    loop_permutation::Vector{Int}
    loop_signs::Vector{Int}
    external_shifts::Vector{Int}
end

function _LoopGCQuotientWorkspace(graph_capacity::Int, loop_capacity::Int)
    return _LoopGCQuotientWorkspace(
        _LoopGCCanonicalizationStorage(graph_capacity),
        _LoopCanonicalGraphBuilder(),
        Vector{Int}(undef, loop_capacity),
        Vector{Int}(undef, loop_capacity),
        Vector{Int}(undef, loop_capacity),
        Vector{Int}(undef, loop_capacity),
        [Tuple{Int,MomentumCoefficient}[] for _ in 1:loop_capacity],
        Vector{Int}(undef, loop_capacity),
        Vector{Int}(undef, loop_capacity),
        Vector{Int}(undef, loop_capacity),
        Vector{Int}(undef, loop_capacity),
    )
end

function _prepare_loop_gc_workspace!(workspace::_LoopGCQuotientWorkspace, nloops::Int)
    empty!(workspace.builder.colors)
    empty!(workspace.builder.edges)

    resize!(workspace.loop_indices, nloops)
    resize!(workspace.pair_vertices, nloops)
    resize!(workspace.positive_vertices, nloops)
    resize!(workspace.negative_vertices, nloops)
    resize!(workspace.ordered_slots, nloops)
    resize!(workspace.loop_permutation, nloops)
    resize!(workspace.loop_signs, nloops)
    resize!(workspace.external_shifts, nloops)

    @inbounds for slot in 1:nloops
        empty!(workspace.loop_incidences[slot])
        workspace.ordered_slots[slot] = slot
        workspace.loop_permutation[slot] = 0
        workspace.loop_signs[slot] = 1
        workspace.external_shifts[slot] = 0
    end
    return workspace
end

function _external_basis_index_noalloc(
    basis::MomentumBasis, external::MomentumVariable
)::Int
    external_index = 0
    @inbounds for basis_index in eachindex(basis.variables)
        basis.variables[basis_index] == external || continue
        iszero(external_index) || throw(
            ArgumentError(
                "external momentum must occur exactly once in the momentum basis"
            ),
        )
        external_index = basis_index
    end
    iszero(external_index) && throw(
        ArgumentError("external momentum must occur exactly once in the momentum basis")
    )
    return external_index
end

function _write_loop_basis_indices!(
    workspace::_LoopGCQuotientWorkspace, basis::MomentumBasis, external::MomentumVariable
)
    external_index = _external_basis_index_noalloc(basis, external)
    loop_slot = 0
    @inbounds for basis_index in eachindex(basis.variables)
        basis_index == external_index && continue
        loop_slot += 1
        workspace.loop_indices[loop_slot] = basis_index
    end
    return external_index
end

@inline function _loop_incidence_vertex_count(
    momentum::LinearMomentum, external_index::Int
)::Int
    count = 0
    @inbounds for basis_index in eachindex(momentum.coefficients)
        basis_index == external_index && continue
        count += !iszero(momentum.coefficients[basis_index])
    end
    return count
end

@inline function _projective_component_vertex_count(component::MomentumComponent)::Int
    nonzero_coefficients = 0
    @inbounds for coefficient in component.momentum.coefficients
        nonzero_coefficients += !iszero(coefficient)
    end
    return 3 + 2 * nonzero_coefficients
end

function _projective_support_vertex_count(
    sector::ReducedCollisionSector, external_index::Int
)::Int
    count = 0
    support = frequency_support(sector)
    for shell in support.shells
        count += 1
        for (atom, _) in shell.energy
            count += 1 + _loop_incidence_vertex_count(atom.momentum, external_index)
        end
    end
    for principal_value in support.principal_values
        count += 1
        for (atom, _) in principal_value.energy
            count += 1 + _loop_incidence_vertex_count(atom.momentum, external_index)
        end
    end
    return count
end

function _loop_gc_capacities(expression::OccupationReducedExpression)
    graph_capacity = 0
    loop_capacity = 0
    for (sector, polynomial) in occupation_reduced_terms(expression)
        basis = momentum_basis(sector)
        external_index = _external_basis_index_noalloc(
            basis, external_wigner_momentum(sector)
        )
        nloops = length(basis) - 1
        loop_capacity = max(loop_capacity, nloops)
        fixed_vertices =
            3 + 3 * nloops + _projective_support_vertex_count(sector, external_index)

        for (occupation_monomial, _) in polynomial
            occupation_vertices = 0
            for atom in occupation_monomial
                occupation_vertices +=
                    1 + _loop_incidence_vertex_count(atom.momentum, external_index)
            end
            for (kinematic_monomial, _) in kinematic_factor(sector)
                kinematic_vertices = 1
                for component in kinematic_monomial
                    kinematic_vertices += _projective_component_vertex_count(component)
                end
                graph_capacity = max(
                    graph_capacity,
                    fixed_vertices + occupation_vertices + kinematic_vertices,
                )
            end
        end
    end
    return graph_capacity, loop_capacity
end

function _order_loop_slots!(
    workspace::_LoopGCQuotientWorkspace,
    buffer::GC.DirectedCanonicalizationBuffer,
    nloops::Int,
)
    @inbounds for index in 2:nloops
        slot = workspace.ordered_slots[index]
        rank = GC.canonical_rank(buffer, workspace.pair_vertices[slot])
        position = index - 1
        while position >= 1
            previous_slot = workspace.ordered_slots[position]
            previous_rank = GC.canonical_rank(
                buffer, workspace.pair_vertices[previous_slot]
            )
            rank < previous_rank || break
            workspace.ordered_slots[position + 1] = previous_slot
            position -= 1
        end
        workspace.ordered_slots[position + 1] = slot
    end
    return workspace.ordered_slots
end

function _graphcombinations_projective_canonical_loop_transform(
    sector::ReducedCollisionSector{S},
    monomial::OccupationMonomial{S},
    workspace::_LoopGCQuotientWorkspace,
) where {S<:Statistics}
    basis = momentum_basis(sector)
    external = external_wigner_momentum(sector)
    nloops = length(basis) - 1
    _prepare_loop_gc_workspace!(workspace, nloops)
    external_index = _write_loop_basis_indices!(workspace, basis, external)

    builder = workspace.builder
    root = _add_loop_vertex!(builder, _loop_graph_color(1))
    for slot in 1:nloops
        pair = _add_loop_vertex!(builder, _loop_graph_color(2))
        positive = _add_loop_vertex!(builder, _loop_graph_color(3))
        negative = _add_loop_vertex!(builder, _loop_graph_color(3))
        workspace.pair_vertices[slot] = pair
        workspace.positive_vertices[slot] = positive
        workspace.negative_vertices[slot] = negative
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
        workspace.loop_indices,
        workspace.positive_vertices,
        workspace.negative_vertices,
        workspace.loop_incidences,
    )

    buffer = _canonicalize_loop_builder!(workspace.canonicalization, builder)
    _order_loop_slots!(workspace, buffer, nloops)
    @inbounds for new_slot in 1:nloops
        old_slot = workspace.ordered_slots[new_slot]
        workspace.loop_permutation[old_slot] = new_slot
        positive_rank = GC.canonical_rank(buffer, workspace.positive_vertices[old_slot])
        negative_rank = GC.canonical_rank(buffer, workspace.negative_vertices[old_slot])
        workspace.loop_signs[old_slot] = positive_rank < negative_rank ? 1 : -1
    end

    return loop_permutation_transform(
        basis,
        external,
        workspace.loop_permutation,
        workspace.loop_signs,
        workspace.external_shifts,
    )
end


function _prepare_matrixfree_gc_gauge!(
    sector::ReducedCollisionSector{S},
    monomial::OccupationMonomial{S},
    workspace::_LoopGCQuotientWorkspace,
) where {S<:Statistics}
    basis = momentum_basis(sector)
    external = external_wigner_momentum(sector)
    nloops = length(basis) - 1
    _prepare_loop_gc_workspace!(workspace, nloops)
    external_index = _write_loop_basis_indices!(workspace, basis, external)

    builder = workspace.builder
    root = _add_loop_vertex!(builder, _loop_graph_color(1))
    for slot in 1:nloops
        pair = _add_loop_vertex!(builder, _loop_graph_color(2))
        positive = _add_loop_vertex!(builder, _loop_graph_color(3))
        negative = _add_loop_vertex!(builder, _loop_graph_color(3))
        workspace.pair_vertices[slot] = pair
        workspace.positive_vertices[slot] = positive
        workspace.negative_vertices[slot] = negative
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
        workspace.loop_indices,
        workspace.positive_vertices,
        workspace.negative_vertices,
        workspace.loop_incidences,
    )

    buffer = _canonicalize_loop_builder!(workspace.canonicalization, builder)
    _order_loop_slots!(workspace, buffer, nloops)
    @inbounds for new_slot in 1:nloops
        old_slot = workspace.ordered_slots[new_slot]
        workspace.loop_permutation[old_slot] = new_slot
        positive_rank = GC.canonical_rank(buffer, workspace.positive_vertices[old_slot])
        negative_rank = GC.canonical_rank(buffer, workspace.negative_vertices[old_slot])
        workspace.loop_signs[old_slot] = positive_rank < negative_rank ? 1 : -1
    end
    return external_index
end

function _matrixfree_gc_momentum(
    momentum::LinearMomentum,
    workspace::_LoopGCQuotientWorkspace,
    external_index::Int,
)
    n = length(momentum)
    n == length(workspace.loop_indices) + 1 || throw(
        DimensionMismatch(
            "momentum and signed loop permutation use different basis sizes"
        ),
    )
    coefficients = Vector{MomentumCoefficient}(undef, n)
    coefficients[external_index] = momentum[external_index]
    @inbounds for old_slot in eachindex(workspace.loop_permutation)
        old_index = workspace.loop_indices[old_slot]
        new_slot = workspace.loop_permutation[old_slot]
        new_index = workspace.loop_indices[new_slot]
        coefficients[new_index] =
            workspace.loop_signs[old_slot] * momentum[old_index]
    end
    return LinearMomentum(coefficients)
end

function _matrixfree_gc_occupation_monomial(
    monomial::OccupationMonomial{S},
    workspace::_LoopGCQuotientWorkspace,
    external_index::Int,
) where {S<:Statistics}
    atoms = OccupationAtom{S}[
        OccupationAtom{S}(
            atom.family,
            _matrixfree_gc_momentum(atom.momentum, workspace, external_index),
        ) for atom in monomial
    ]
    return OccupationMonomial(atoms)
end

function _matrixfree_gc_momentum_monomial(
    monomial::MomentumMonomial,
    workspace::_LoopGCQuotientWorkspace,
    external_index::Int,
)
    factors = MomentumComponent[
        MomentumComponent(
            _matrixfree_gc_momentum(
                component.momentum, workspace, external_index
            ),
            component.axis,
        ) for component in monomial
    ]
    return MomentumMonomial(factors)
end

function _matrixfree_gc_momentum_polynomial(
    polynomial::MomentumPolynomial{C},
    workspace::_LoopGCQuotientWorkspace,
    external_index::Int,
) where {C<:Number}
    terms = Pair{MomentumMonomial,C}[]
    sizehint!(terms, length(polynomial))
    for (monomial, coefficient) in polynomial
        transformed = _matrixfree_gc_momentum_monomial(
            monomial, workspace, external_index
        )
        normalized, factor, nonzero = _projective_kinematic_monomial(transformed)
        nonzero || continue
        push!(terms, normalized => coefficient * convert(C, factor))
    end
    return MomentumPolynomial{C}(terms)
end

function _matrixfree_gc_energy_form(
    form::EnergyForm{S},
    workspace::_LoopGCQuotientWorkspace,
    external_index::Int,
) where {S<:Statistics}
    energy_basis_size(form) == length(workspace.loop_indices) + 1 || throw(
        DimensionMismatch(
            "energy form and signed loop permutation use different basis sizes"
        ),
    )
    terms = Pair{DispersionAtom{S},EnergyCoefficient}[
        DispersionAtom{S}(
            atom.family,
            _matrixfree_gc_momentum(atom.momentum, workspace, external_index),
        ) => coefficient for (atom, coefficient) in form
    ]
    return EnergyForm(energy_basis_size(form), terms)
end

function _matrixfree_gc_frequency_support(
    support::FrequencySupport{S},
    workspace::_LoopGCQuotientWorkspace,
    external_index::Int,
) where {S<:Statistics}
    shells = EnergyShell{S}[]
    principal_values = PrincipalValueSupport{S}[]
    factor = one(MomentumCoefficient)
    sizehint!(shells, length(support.shells))
    sizehint!(principal_values, length(support.principal_values))

    for shell in support.shells
        transformed, shell_factor = energy_shell(
            _matrixfree_gc_energy_form(shell.energy, workspace, external_index)
        )
        push!(shells, transformed)
        factor *= shell_factor
    end
    for principal_value in support.principal_values
        transformed, pv_factor = principal_value_support(
            _matrixfree_gc_energy_form(
                principal_value.energy, workspace, external_index
            )
        )
        push!(principal_values, transformed)
        factor *= pv_factor
    end
    return FrequencySupport(shells, principal_values), factor
end

function _matrixfree_gc_sector(
    sector::ReducedCollisionSector{S},
    workspace::_LoopGCQuotientWorkspace,
    external_index::Int,
) where {S<:Statistics}
    support, factor = _matrixfree_gc_frequency_support(
        frequency_support(sector), workspace, external_index
    )
    transformed = CollisionKernelSector{S}(
        parameters(sector),
        momentum_basis(sector),
        external_wigner_momentum(sector),
        _matrixfree_gc_momentum_polynomial(
            kinematic_factor(sector), workspace, external_index
        ),
        support,
    )
    return transformed, factor
end

function _quotient_loop_momenta_graphcombinations(
    expression::OccupationReducedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:Statistics,O,G,Ctx<:AbstractWignerContext}
    D = promote_type(C, ComplexRationals, Rational{Int})
    out = Dict{CollisionKernelSector{S},OccupationPolynomial{D,S}}()
    graph_capacity, loop_capacity = _loop_gc_capacities(expression)
    workspace = _LoopGCQuotientWorkspace(graph_capacity, loop_capacity)

    for (sector, polynomial) in occupation_reduced_terms(expression)
        for (occupation_monomial, occupation_coefficient) in polynomial
            for (kinematic_monomial, kinematic_coefficient) in kinematic_factor(sector)
                atom_sector = _kinematic_atom_sector(sector, kinematic_monomial)
                external_index = _prepare_matrixfree_gc_gauge!(
                    atom_sector, occupation_monomial, workspace
                )
                transformed_sector, support_factor = _matrixfree_gc_sector(
                    atom_sector, workspace, external_index
                )
                transformed_monomial = _matrixfree_gc_occupation_monomial(
                    occupation_monomial, workspace, external_index
                )
                transformed_coefficient =
                    convert(D, occupation_coefficient) *
                    convert(D, kinematic_coefficient) *
                    convert(D, support_factor)
                contribution = Pair{OccupationMonomial{S},D}[
                    transformed_monomial => transformed_coefficient
                ]
                _push_kernel_polynomial!(
                    out,
                    transformed_sector,
                    OccupationPolynomial{D,S}(contribution),
                )
            end
        end
    end

    return LoopQuotientedExpression{D,S,O,G,Ctx}(
        out,
        target_family(expression),
        parameters(expression),
        wigner_context(expression),
    )
end

function quotient_loop_momenta(
    expression::OccupationReducedExpression{C,Boson,O,G,Ctx}
) where {C<:Number,O,G,Ctx<:AbstractWignerContext}
    return _quotient_loop_momenta_graphcombinations(expression)
end

function quotient_loop_momenta(
    expression::OccupationReducedExpression{C,Fermion,O,G,Ctx}
) where {C<:Number,O,G,Ctx<:AbstractWignerContext}
    return _quotient_loop_momenta_graphcombinations(expression)
end
