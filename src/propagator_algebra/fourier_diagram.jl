"""
Momentum-space identity for a coordinate-space `Diagram`.

The coordinate diagram remains the canonical graph/field identity. `edge_momenta`
contains the exact momentum assigned to each coordinate-space edge in the same fixed
edge order. The first `external_count` variables in `basis` are external momenta and
the remaining `loop_count` variables are canonical loop momenta. `kinematic` is kept
separate from the numerical diagram coefficient; this routing-only layer uses
`nothing` until derivative momentum polynomials are introduced.
"""
struct FourierDiagram{S<:Statistics,E1,E2,K}
    coordinate::Diagram{S,E1,E2}
    basis::MomentumBasis
    edge_momenta::FixedVector{E1,LinearMomentum}
    external_count::Int16
    loop_count::Int16
    kinematic::K
end

coordinate_diagram(diagram::FourierDiagram) = diagram.coordinate
momentum_basis(diagram::FourierDiagram) = diagram.basis
edge_momenta(diagram::FourierDiagram) = diagram.edge_momenta
external_momentum_count(diagram::FourierDiagram) = Int(diagram.external_count)
loop_momentum_count(diagram::FourierDiagram) = Int(diagram.loop_count)
kinematic_factor(diagram::FourierDiagram) = diagram.kinematic

function Base.isequal(
    a::FourierDiagram{S,E1,E2,K1}, b::FourierDiagram{S,E1,E2,K2}
) where {S,E1,E2,K1,K2}
    return isequal(a.coordinate, b.coordinate) &&
           isequal(a.basis, b.basis) &&
           isequal(a.edge_momenta, b.edge_momenta) &&
           a.external_count == b.external_count &&
           a.loop_count == b.loop_count &&
           isequal(a.kinematic, b.kinematic)
end
Base.:(==)(a::FourierDiagram, b::FourierDiagram) = isequal(a, b)
function Base.hash(x::FourierDiagram, h::UInt)
    h = hash(FourierDiagram, h)
    h = hash(x.coordinate, h)
    h = hash(x.basis, h)
    h = hash(x.edge_momenta, h)
    h = hash(x.external_count, h)
    h = hash(x.loop_count, h)
    return hash(x.kinematic, h)
end

@inline function _routing_edge_endpoints(c::Contraction)
    out_position, in_position = positions(c)
    return convert(Int16, index(in_position)), convert(Int16, index(out_position))
end

function _canonical_internal_routing_edges(legged::Vector{Contraction{S}}) where {S<:Statistics}
    counts = Dict{Tuple{Int16,Int16},Int}()
    routing_edges = RoutingEdge[]
    legged_indices = Int[]

    for (i, contraction) in enumerate(legged)
        out_position, in_position = positions(contraction)
        if is_bulk(out_position) && is_bulk(in_position)
            tail, head = _routing_edge_endpoints(contraction)
            key = (tail, head)
            ordinal = get(counts, key, 0) + 1
            counts[key] = ordinal
            push!(routing_edges, RoutingEdge(tail, head, ordinal))
            push!(legged_indices, i)
        end
    end

    permutation = sortperm(routing_edges)
    return routing_edges[permutation], legged_indices[permutation]
end

function _bulk_vertices(legged::Vector{Contraction{S}}) where {S<:Statistics}
    vertices = Int16[]
    for contraction in legged
        for position in positions(contraction)
            is_bulk(position) && push!(vertices, convert(Int16, index(position)))
        end
    end
    sort!(unique!(vertices))
    isempty(vertices) && throw(ArgumentError("Fourier routing requires at least one bulk vertex"))
    return vertices
end

function _validate_two_point_external_flow(legged::Vector{Contraction{S}}) where {S<:Statistics}
    incoming = 0
    outgoing = 0
    for contraction in legged
        out_position, in_position = positions(contraction)
        has_in_edge = is_in(out_position) || is_in(in_position)
        has_out_edge = is_out(out_position) || is_out(in_position)
        if has_in_edge || has_out_edge
            (has_in_edge ⊻ has_out_edge) || throw(
                ArgumentError("each external propagator must connect one bulk and one external endpoint")
            )
            count(is_bulk, (out_position, in_position)) == 1 || throw(
                ArgumentError("each external propagator must connect exactly one bulk endpoint")
            )
        end
        if has_in_edge
            is_in(in_position) || throw(
                ArgumentError("incoming external momentum must enter through an edge in-endpoint")
            )
            incoming += 1
        end
        if has_out_edge
            is_out(out_position) || throw(
                ArgumentError("outgoing external momentum must leave through an edge out-endpoint")
            )
            outgoing += 1
        end
    end
    incoming == 1 || throw(ArgumentError("expected exactly one incoming external leg"))
    outgoing == 1 || throw(ArgumentError("expected exactly one outgoing external leg"))
    return nothing
end

function _diagram_affine_routing(d::Diagram{S,E1,E2}) where {S<:Statistics,E1,E2}
    legged, added = with_external_legs(d)
    added >= 0 || throw(ArgumentError("cannot infer two-point external attachments"))
    _validate_two_point_external_flow(legged)

    vertices = _bulk_vertices(legged)
    vertex_index = Dict{Int16,Int}(vertex => i for (i, vertex) in enumerate(vertices))
    routing_edges, legged_indices = _canonical_internal_routing_edges(legged)

    incidence = zeros(Int, length(vertices), length(routing_edges))
    for (column, edge) in enumerate(routing_edges)
        incidence[vertex_index[edge.tail], column] -= 1
        incidence[vertex_index[edge.head], column] += 1
    end

    source = zeros(Int, length(vertices), 1)
    for contraction in legged
        out_position, in_position = positions(contraction)
        has_external =
            is_in(out_position) || is_out(out_position) || is_in(in_position) || is_out(in_position)
        has_external || continue

        if is_bulk(out_position)
            source[vertex_index[convert(Int16, index(out_position))], 1] -= 1
        elseif is_bulk(in_position)
            source[vertex_index[convert(Int16, index(in_position))], 1] += 1
        else
            throw(ArgumentError("external propagator has no bulk endpoint"))
        end
    end

    routing = exact_affine_momentum_routing(incidence, source)
    external_count = external_momentum_count(routing)
    external_count == 1 || error("two-point Fourier routing must have one external momentum")
    loop_count = loop_momentum_count(routing)

    assigned = Vector{LinearMomentum}(undef, E1)
    external_momentum = basis_momentum(routing.basis, 1)

    for i in 1:E1
        out_position, in_position = positions(legged[i])
        if is_bulk(out_position) && is_bulk(in_position)
            continue
        end
        assigned[i] = external_momentum
    end

    for (column, legged_index) in enumerate(legged_indices)
        legged_index <= E1 || error("virtual external-leg reconstruction produced an internal edge")
        assigned[legged_index] = routing.edge_momenta[column]
    end

    fixed_momenta = FixedVector{E1,LinearMomentum}(assigned)
    return routing.basis,
    fixed_momenta,
    convert(Int16, external_count),
    convert(Int16, loop_count)
end

function FourierDiagram(d::Diagram{S,E1,E2}) where {S<:Statistics,E1,E2}
    basis, momenta, external_count, loop_count = _diagram_affine_routing(d)
    return FourierDiagram{S,E1,E2,Nothing}(
        d, basis, momenta, external_count, loop_count, nothing
    )
end
