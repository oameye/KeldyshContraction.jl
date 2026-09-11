"""
Canonical directed graph edge used by the routing layer.

`ordinal` distinguishes parallel edges with the same endpoints. It must be assigned
from a deterministic upstream edge identity; zero is the convenience value for a
simple edge when no parallel-edge distinction is needed.
"""
struct RoutingEdge
    tail::Int16
    head::Int16
    ordinal::Int16
end

function RoutingEdge(tail::Integer, head::Integer, ordinal::Integer=0)
    typemin(Int16) <= tail <= typemax(Int16) ||
        throw(ArgumentError("routing-edge tail is outside Int16 range"))
    typemin(Int16) <= head <= typemax(Int16) ||
        throw(ArgumentError("routing-edge head is outside Int16 range"))
    0 <= ordinal <= typemax(Int16) ||
        throw(ArgumentError("routing-edge ordinal must be non-negative and fit Int16"))
    return RoutingEdge(convert(Int16, tail), convert(Int16, head), convert(Int16, ordinal))
end

function Base.isless(a::RoutingEdge, b::RoutingEdge)
    a.tail == b.tail || return a.tail < b.tail
    a.head == b.head || return a.head < b.head
    return a.ordinal < b.ordinal
end
function Base.isequal(a::RoutingEdge, b::RoutingEdge)
    return a.tail == b.tail && a.head == b.head && a.ordinal == b.ordinal
end
Base.:(==)(a::RoutingEdge, b::RoutingEdge) = isequal(a, b)
function Base.hash(x::RoutingEdge, h::UInt)
    return hash(RoutingEdge, hash(x.ordinal, hash(x.head, hash(x.tail, h))))
end

"""
Exact solution of a momentum-conservation system.

`edge_momenta[i]` is the linear momentum on conservation-matrix column `i`.
`pivot_columns` and `free_columns` refer to that canonical column ordering. For
routing constructed from a graph, `edges` contains the corresponding canonical
edge labels; for a bare conservation matrix it is empty because no graph labels
were supplied.
"""
struct MomentumRouting
    edges::Vector{RoutingEdge}
    basis::MomentumBasis
    edge_momenta::Vector{LinearMomentum}
    pivot_columns::Vector{Int}
    free_columns::Vector{Int}
end

Base.length(routing::MomentumRouting) = length(routing.edge_momenta)
function Base.isequal(a::MomentumRouting, b::MomentumRouting)
    return isequal(a.edges, b.edges) &&
           isequal(a.basis, b.basis) &&
           isequal(a.edge_momenta, b.edge_momenta) &&
           isequal(a.pivot_columns, b.pivot_columns) &&
           isequal(a.free_columns, b.free_columns)
end
Base.:(==)(a::MomentumRouting, b::MomentumRouting) = isequal(a, b)
function Base.hash(x::MomentumRouting, h::UInt)
    h = hash(MomentumRouting, h)
    h = hash(x.edges, h)
    h = hash(x.basis, h)
    h = hash(x.edge_momenta, h)
    h = hash(x.pivot_columns, h)
    return hash(x.free_columns, h)
end

"""Return the edge-by-basis coefficient matrix of an exact routing."""
function routing_matrix(routing::MomentumRouting)::Matrix{MomentumCoefficient}
    nedge = length(routing.edge_momenta)
    nbasis = length(routing.basis)
    out = zeros(MomentumCoefficient, nedge, nbasis)
    for edge in 1:nedge
        coefficients = routing.edge_momenta[edge].coefficients
        length(coefficients) == nbasis || error("inconsistent momentum-routing basis")
        for basis in 1:nbasis
            out[edge, basis] = coefficients[basis]
        end
    end
    return out
end

function _exact_rref(A::Matrix{MomentumCoefficient})
    M = copy(A)
    nrows, ncols = size(M)
    pivots = Int[]
    row = 1

    for col in 1:ncols
        row > nrows && break

        pivot = 0
        for candidate in row:nrows
            if !iszero(M[candidate, col])
                pivot = candidate
                break
            end
        end
        iszero(pivot) && continue

        if pivot != row
            for j in 1:ncols
                M[row, j], M[pivot, j] = M[pivot, j], M[row, j]
            end
        end

        pivot_value = M[row, col]
        for j in col:ncols
            M[row, j] /= pivot_value
        end

        for other in 1:nrows
            other == row && continue
            factor = M[other, col]
            iszero(factor) && continue
            for j in col:ncols
                M[other, j] -= factor * M[row, j]
            end
        end

        push!(pivots, col)
        row += 1
    end

    return M, pivots
end

"""
Construct the canonical exact nullspace routing of an integer conservation matrix.

Columns are edge momenta and rows are conservation equations. Pivot columns are
chosen left-to-right. Free columns, also left-to-right, define independent momentum
variables with coefficient `+1`; all dependent edge momenta are then exact rational
linear combinations of those variables.
"""
function exact_momentum_routing(A::AbstractMatrix{<:Integer})::MomentumRouting
    nrows, ncols = size(A)
    exact = Matrix{MomentumCoefficient}(undef, nrows, ncols)
    for j in 1:ncols, i in 1:nrows
        exact[i, j] = convert(Int, A[i, j]) // 1
    end

    reduced, pivots = _exact_rref(exact)
    pivot_mask = falses(ncols)
    for col in pivots
        pivot_mask[col] = true
    end
    free = Int[col for col in 1:ncols if !pivot_mask[col]]
    basis = MomentumBasis(length(free))

    free_to_basis = zeros(Int, ncols)
    for (basis_index, col) in enumerate(free)
        free_to_basis[col] = basis_index
    end

    pivot_to_row = zeros(Int, ncols)
    for (pivot_row, col) in enumerate(pivots)
        pivot_to_row[col] = pivot_row
    end

    momenta = Vector{LinearMomentum}(undef, ncols)
    for col in 1:ncols
        coefficients = zeros(MomentumCoefficient, length(free))
        basis_index = free_to_basis[col]
        if !iszero(basis_index)
            coefficients[basis_index] = 1 // 1
        else
            pivot_row = pivot_to_row[col]
            for (j, free_col) in enumerate(free)
                coefficients[j] = -reduced[pivot_row, free_col]
            end
        end
        momenta[col] = LinearMomentum(coefficients)
    end

    return MomentumRouting(RoutingEdge[], basis, momenta, pivots, free)
end

function _canonical_incidence(edges::Vector{RoutingEdge})
    canonical_edges = sort(copy(edges))
    for i in 2:length(canonical_edges)
        isequal(canonical_edges[i - 1], canonical_edges[i]) && throw(
            ArgumentError(
                "parallel routing edges must have distinct deterministic ordinals"
            ),
        )
    end

    vertices = sort!(
        unique(
            Int16[vertex for edge in canonical_edges for vertex in (edge.tail, edge.head)]
        ),
    )
    vertex_index = Dict{Int16,Int}(vertex => i for (i, vertex) in enumerate(vertices))
    incidence = zeros(Int, length(vertices), length(canonical_edges))

    for (column, edge) in enumerate(canonical_edges)
        tail = vertex_index[edge.tail]
        head = vertex_index[edge.head]
        incidence[tail, column] -= 1
        incidence[head, column] += 1
    end
    return canonical_edges, incidence
end

"""
Route momenta on a directed graph using exact incidence-matrix nullspace algebra.

Input edge order is irrelevant: edges are canonicalized lexicographically by
`(tail, head, ordinal)` before constructing the incidence matrix. Parallel edges
must carry distinct deterministic ordinals so the routing remains attachable to a
later Fourier-diagram representation. Self-loops become zero incidence columns and
correctly contribute independent cycle variables.
"""
function momentum_routing(edges::Vector{RoutingEdge})::MomentumRouting
    canonical_edges, incidence = _canonical_incidence(edges)
    raw = exact_momentum_routing(incidence)
    return MomentumRouting(
        canonical_edges, raw.basis, raw.edge_momenta, raw.pivot_columns, raw.free_columns
    )
end
