# Production GraphCombinations backend for the physical dummy-loop quotient.
#
# The Nauty projective transform in `collision_momentum_projective.jl` remains an
# independent private oracle. Public Boson/Fermion quotienting is specialized below
# to use the exact GraphCombinations directed-canonicalization workspace.

struct _LoopGCCanonicalizationBuffers
    workspace::GraphCombinations.DirectedCanonicalizationWorkspace
    buffer::GraphCombinations.DirectedCanonicalizationBuffer
end

const _LoopGCCanonicalizationCache = Dict{Int,_LoopGCCanonicalizationBuffers}

function _loop_gc_graph(builder::_LoopCanonicalGraphBuilder)
    color_classes = sort!(unique(copy(builder.colors)))
    labels = Int[searchsortedfirst(color_classes, color) for color in builder.colors]
    edges = Pair{Int,Int}[source => target for (source, target) in builder.edges]
    return GraphCombinations.DirectedGCGraph(edges, length(labels)), labels
end

function _loop_gc_buffers!(cache::_LoopGCCanonicalizationCache, num_vertices::Int)
    if haskey(cache, num_vertices)
        return cache[num_vertices]
    end
    buffers = _LoopGCCanonicalizationBuffers(
        GraphCombinations.DirectedCanonicalizationWorkspace(num_vertices),
        GraphCombinations.DirectedCanonicalizationBuffer(num_vertices),
    )
    cache[num_vertices] = buffers
    return buffers
end

function _graphcombinations_projective_canonical_loop_transform(
    sector::ReducedCollisionSector{S},
    monomial::OccupationMonomial{S},
    cache::_LoopGCCanonicalizationCache,
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

    graph, labels = _loop_gc_graph(builder)
    buffers = _loop_gc_buffers!(cache, graph.num_vertices)
    GraphCombinations.canonicalize_directed!(
        buffers.buffer, buffers.workspace, graph, labels
    )

    ordered_slots = sortperm(
        1:nloops;
        by=slot -> GraphCombinations.canonical_rank(buffers.buffer, pair_vertices[slot]),
    )
    loop_permutation = zeros(Int, nloops)
    loop_signs = ones(Int, nloops)
    for (new_slot, old_slot) in enumerate(ordered_slots)
        loop_permutation[old_slot] = new_slot
        positive_rank = GraphCombinations.canonical_rank(
            buffers.buffer, positive_vertices[old_slot]
        )
        negative_rank = GraphCombinations.canonical_rank(
            buffers.buffer, negative_vertices[old_slot]
        )
        loop_signs[old_slot] = positive_rank < negative_rank ? 1 : -1
    end

    return loop_permutation_transform(
        basis, external, loop_permutation, loop_signs, zeros(Int, nloops)
    )
end

function _quotient_loop_momenta_graphcombinations(
    expression::OccupationReducedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:Statistics,O,G,Ctx<:AbstractWignerContext}
    D = promote_type(C, ComplexRationals, Rational{Int})
    out = Dict{CollisionKernelSector{S},OccupationPolynomial{D,S}}()
    cache = _LoopGCCanonicalizationCache()

    for (sector, polynomial) in occupation_reduced_terms(expression)
        for (occupation_monomial, occupation_coefficient) in polynomial
            for (kinematic_monomial, kinematic_coefficient) in kinematic_factor(sector)
                atom_sector = _kinematic_atom_sector(sector, kinematic_monomial)
                transform = _graphcombinations_projective_canonical_loop_transform(
                    atom_sector, occupation_monomial, cache
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
