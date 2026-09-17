# Production GraphCombinations backend for the physical dummy-loop quotient.
#
# The Nauty projective transform in `collision_momentum_projective.jl` remains an
# independent private oracle. Public Boson/Fermion quotienting is specialized below
# to use the exact GraphCombinations directed-canonicalization workspace.

struct _LoopGCCanonicalizationBuffers
    workspace::GC.DirectedCanonicalizationWorkspace
    buffer::GC.DirectedCanonicalizationBuffer
    graph::GC.DirectedGCGraph
    labels::Vector{Int}
    color_order::Vector{_LoopGraphColor}
end

const _LoopGCCanonicalizationCache = Dict{Int,_LoopGCCanonicalizationBuffers}

mutable struct _LoopGCQuotientWorkspace
    cache::_LoopGCCanonicalizationCache
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

function _LoopGCQuotientWorkspace()
    return _LoopGCQuotientWorkspace(
        _LoopGCCanonicalizationCache(),
        _LoopCanonicalGraphBuilder(),
        Int[],
        Int[],
        Int[],
        Int[],
        Vector{Vector{Tuple{Int,MomentumCoefficient}}}(),
        Int[],
        Int[],
        Int[],
        Int[],
    )
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

function _prepare_loop_gc_graph!(
    buffers::_LoopGCCanonicalizationBuffers, builder::_LoopCanonicalGraphBuilder
)
    num_vertices = length(builder.colors)
    buffers.graph.num_vertices == num_vertices ||
        throw(DimensionMismatch("cached directed graph has the wrong size"))

    copyto!(buffers.color_order, builder.colors)
    sort!(buffers.color_order)
    @inbounds for vertex in eachindex(builder.colors)
        buffers.labels[vertex] = searchsortedfirst(
            buffers.color_order, builder.colors[vertex]
        )
    end

    fill!(buffers.graph.multiplicities, 0)
    @inbounds for (source, target) in builder.edges
        slot = (source - 1) * num_vertices + target
        buffers.graph.multiplicities[slot] += 1
    end
    return buffers
end

function _loop_gc_buffers!(cache::_LoopGCCanonicalizationCache, num_vertices::Int)
    if haskey(cache, num_vertices)
        return cache[num_vertices]
    end
    buffers = _LoopGCCanonicalizationBuffers(num_vertices)
    cache[num_vertices] = buffers
    return buffers
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

    while length(workspace.loop_incidences) < nloops
        push!(workspace.loop_incidences, Tuple{Int,MomentumCoefficient}[])
    end
    @inbounds for slot in 1:nloops
        empty!(workspace.loop_incidences[slot])
        workspace.ordered_slots[slot] = slot
        workspace.loop_permutation[slot] = 0
        workspace.loop_signs[slot] = 1
        workspace.external_shifts[slot] = 0
    end
    return workspace
end

function _write_loop_basis_indices!(
    workspace::_LoopGCQuotientWorkspace,
    basis::MomentumBasis,
    external::MomentumVariable,
)
    external_index = 0
    @inbounds for basis_index in eachindex(basis.variables)
        basis.variables[basis_index] == external || continue
        iszero(external_index) || throw(
            ArgumentError("external momentum must occur exactly once in the momentum basis")
        )
        external_index = basis_index
    end
    iszero(external_index) && throw(
        ArgumentError("external momentum must occur exactly once in the momentum basis")
    )

    loop_slot = 0
    @inbounds for basis_index in eachindex(basis.variables)
        basis_index == external_index && continue
        loop_slot += 1
        workspace.loop_indices[loop_slot] = basis_index
    end
    return external_index
end

function _order_loop_slots!(
    workspace::_LoopGCQuotientWorkspace,
    buffers::_LoopGCCanonicalizationBuffers,
    nloops::Int,
)
    @inbounds for index in 2:nloops
        slot = workspace.ordered_slots[index]
        rank = GC.canonical_rank(buffers.buffer, workspace.pair_vertices[slot])
        position = index - 1
        while position >= 1
            previous_slot = workspace.ordered_slots[position]
            previous_rank = GC.canonical_rank(
                buffers.buffer, workspace.pair_vertices[previous_slot]
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

    buffers = _loop_gc_buffers!(workspace.cache, length(builder.colors))
    _prepare_loop_gc_graph!(buffers, builder)
    GC.canonicalize_directed!(
        buffers.buffer, buffers.workspace, buffers.graph, buffers.labels
    )

    _order_loop_slots!(workspace, buffers, nloops)
    @inbounds for new_slot in 1:nloops
        old_slot = workspace.ordered_slots[new_slot]
        workspace.loop_permutation[old_slot] = new_slot
        positive_rank = GC.canonical_rank(
            buffers.buffer, workspace.positive_vertices[old_slot]
        )
        negative_rank = GC.canonical_rank(
            buffers.buffer, workspace.negative_vertices[old_slot]
        )
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

function _quotient_loop_momenta_graphcombinations(
    expression::OccupationReducedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:Statistics,O,G,Ctx<:AbstractWignerContext}
    D = promote_type(C, ComplexRationals, Rational{Int})
    out = Dict{CollisionKernelSector{S},OccupationPolynomial{D,S}}()
    workspace = _LoopGCQuotientWorkspace()

    for (sector, polynomial) in occupation_reduced_terms(expression)
        for (occupation_monomial, occupation_coefficient) in polynomial
            for (kinematic_monomial, kinematic_coefficient) in kinematic_factor(sector)
                atom_sector = _kinematic_atom_sector(sector, kinematic_monomial)
                transform = _graphcombinations_projective_canonical_loop_transform(
                    atom_sector, occupation_monomial, workspace
                )
                transformed_sector, support_factor = _transform_kernel_sector(
                    atom_sector, transform
                )
                transformed_monomial = transform_loop_momenta(
                    occupation_monomial, transform
                )
                transformed_coefficient =
                    convert(D, occupation_coefficient) *
                    convert(D, kinematic_coefficient) *
                    convert(D, support_factor)
                contribution = Pair{OccupationMonomial{S},D}[
                    transformed_monomial => transformed_coefficient
                ]
                _push_kernel_polynomial!(
                    out, transformed_sector, OccupationPolynomial{D,S}(contribution)
                )
            end
        end
    end

    return LoopQuotientedExpression{D,S,O,G,Ctx}(
        out, target_family(expression), parameters(expression), wigner_context(expression)
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
