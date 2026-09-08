"""Exact affine momentum routing with external and loop basis variables."""
struct AffineMomentumRouting
    basis::MomentumBasis
    edge_momenta::Vector{LinearMomentum}
    external_count::Int16
    pivot_columns::Vector{Int}
    free_columns::Vector{Int}
end

Base.length(routing::AffineMomentumRouting) = length(routing.edge_momenta)

function Base.isequal(a::AffineMomentumRouting, b::AffineMomentumRouting)
    return isequal(a.basis, b.basis) &&
           isequal(a.edge_momenta, b.edge_momenta) &&
           a.external_count == b.external_count &&
           isequal(a.pivot_columns, b.pivot_columns) &&
           isequal(a.free_columns, b.free_columns)
end
Base.:(==)(a::AffineMomentumRouting, b::AffineMomentumRouting) = isequal(a, b)
function Base.hash(x::AffineMomentumRouting, h::UInt)
    h = hash(AffineMomentumRouting, h)
    h = hash(x.basis, h)
    h = hash(x.edge_momenta, h)
    h = hash(x.external_count, h)
    h = hash(x.pivot_columns, h)
    return hash(x.free_columns, h)
end

external_momentum_count(routing::AffineMomentumRouting) = Int(routing.external_count)
function loop_momentum_count(routing::AffineMomentumRouting)
    return length(routing.basis) - external_momentum_count(routing)
end

"""Return the edge-by-basis coefficient matrix of an affine routing."""
function routing_matrix(routing::AffineMomentumRouting)::Matrix{MomentumCoefficient}
    nedge = length(routing.edge_momenta)
    nbasis = length(routing.basis)
    out = zeros(MomentumCoefficient, nedge, nbasis)
    for edge in 1:nedge
        coefficients = routing.edge_momenta[edge].coefficients
        length(coefficients) == nbasis ||
            error("inconsistent affine momentum-routing basis")
        for basis in 1:nbasis
            out[edge, basis] = coefficients[basis]
        end
    end
    return out
end

function _exact_integer_matrix(A::AbstractMatrix{<:Integer})::Matrix{MomentumCoefficient}
    nrows, ncols = size(A)
    exact = Matrix{MomentumCoefficient}(undef, nrows, ncols)
    for j in 1:ncols, i in 1:nrows
        exact[i, j] = convert(Int, A[i, j]) // 1
    end
    return exact
end

function _exact_rref_with_source(
    A::Matrix{MomentumCoefficient}, source::Matrix{MomentumCoefficient}
)
    nrows, ncols = size(A)
    size(source, 1) == nrows || throw(DimensionMismatch("source row count must match A"))
    nsource = size(source, 2)
    total = ncols + nsource
    augmented = Matrix{MomentumCoefficient}(undef, nrows, total)

    for j in 1:ncols, i in 1:nrows
        augmented[i, j] = A[i, j]
    end
    for j in 1:nsource, i in 1:nrows
        augmented[i, ncols + j] = source[i, j]
    end

    pivots = Int[]
    row = 1
    for col in 1:ncols
        row > nrows && break

        pivot = 0
        for candidate in row:nrows
            if !iszero(augmented[candidate, col])
                pivot = candidate
                break
            end
        end
        iszero(pivot) && continue

        if pivot != row
            for j in 1:total
                augmented[row, j], augmented[pivot, j] = augmented[pivot, j],
                augmented[row, j]
            end
        end

        pivot_value = augmented[row, col]
        for j in col:total
            augmented[row, j] /= pivot_value
        end

        for other in 1:nrows
            other == row && continue
            factor = augmented[other, col]
            iszero(factor) && continue
            for j in col:total
                augmented[other, j] -= factor * augmented[row, j]
            end
        end

        push!(pivots, col)
        row += 1
    end

    reduced = Matrix{MomentumCoefficient}(undef, nrows, ncols)
    transformed_source = Matrix{MomentumCoefficient}(undef, nrows, nsource)
    for j in 1:ncols, i in 1:nrows
        reduced[i, j] = augmented[i, j]
    end
    for j in 1:nsource, i in 1:nrows
        transformed_source[i, j] = augmented[i, ncols + j]
    end
    return reduced, transformed_source, pivots
end

function _check_affine_consistency(
    reduced::Matrix{MomentumCoefficient}, source::Matrix{MomentumCoefficient}
)::Nothing
    nrows, ncols = size(reduced)
    nsource = size(source, 2)
    for row in 1:nrows
        has_lhs = false
        for col in 1:ncols
            if !iszero(reduced[row, col])
                has_lhs = true
                break
            end
        end
        has_lhs && continue
        for col in 1:nsource
            iszero(source[row, col]) ||
                throw(ArgumentError("inconsistent affine momentum-conservation system"))
        end
    end
    return nothing
end

"""
Construct an exact affine solution of `A * p = source * q`.

The first `size(source, 2)` variables in the returned basis are external momenta.
Free columns of `A` then define loop momenta, normalized left-to-right with coefficient
`+1`. All coefficients are exact `Rational{Int}` values.
"""
function exact_affine_momentum_routing(
    A::AbstractMatrix{<:Integer}, source::AbstractMatrix{<:Integer}
)::AffineMomentumRouting
    size(A, 1) == size(source, 1) ||
        throw(DimensionMismatch("source row count must match conservation matrix"))

    exact = _exact_integer_matrix(A)
    exact_source = _exact_integer_matrix(source)
    reduced, transformed_source, pivots = _exact_rref_with_source(exact, exact_source)
    _check_affine_consistency(reduced, transformed_source)

    ncols = size(A, 2)
    nsource = size(source, 2)
    nsource <= typemax(Int16) ||
        throw(ArgumentError("too many external momentum variables"))

    pivot_mask = falses(ncols)
    for col in pivots
        pivot_mask[col] = true
    end
    free = Int[col for col in 1:ncols if !pivot_mask[col]]
    basis = MomentumBasis(nsource + length(free))

    free_to_basis = zeros(Int, ncols)
    for (loop_index, col) in enumerate(free)
        free_to_basis[col] = nsource + loop_index
    end

    pivot_to_row = zeros(Int, ncols)
    for (pivot_row, col) in enumerate(pivots)
        pivot_to_row[col] = pivot_row
    end

    momenta = Vector{LinearMomentum}(undef, ncols)
    for col in 1:ncols
        coefficients = zeros(MomentumCoefficient, length(basis))
        basis_index = free_to_basis[col]
        if !iszero(basis_index)
            coefficients[basis_index] = 1 // 1
        else
            pivot_row = pivot_to_row[col]
            for external in 1:nsource
                coefficients[external] = transformed_source[pivot_row, external]
            end
            for (loop_index, free_col) in enumerate(free)
                coefficients[nsource + loop_index] = -reduced[pivot_row, free_col]
            end
        end
        momenta[col] = LinearMomentum(coefficients)
    end

    return AffineMomentumRouting(basis, momenta, convert(Int16, nsource), pivots, free)
end
