using KeldyshContraction, Test
using KeldyshContraction:
    LinearMomentum,
    MomentumBasis,
    MomentumRouting,
    MomentumVariable,
    RoutingEdge,
    basis_momentum,
    exact_momentum_routing,
    momentum_routing,
    routing_matrix

function routing_recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        routing_recursively_concrete(eltype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        routing_recursively_concrete(FT, seen) || return false
    end
    return true
end

function connected_cycle_rank(edges::Vector{RoutingEdge})
    vertices = unique(Int16[v for edge in edges for v in (edge.tail, edge.head)])
    return length(edges) - length(vertices) + 1
end

function independent_incidence(edges::Vector{RoutingEdge})
    canonical = sort(copy(edges))
    vertices = sort!(unique(Int16[v for edge in canonical for v in (edge.tail, edge.head)]))
    index = Dict(v => i for (i, v) in enumerate(vertices))
    incidence = zeros(Int, length(vertices), length(canonical))
    for (column, edge) in enumerate(canonical)
        incidence[index[edge.tail], column] -= 1
        incidence[index[edge.head], column] += 1
    end
    return incidence
end

@testset "exact momentum algebra" begin
    basis = @inferred MomentumBasis(3)
    @test length(basis) == 3
    @test basis[2] == MomentumVariable(2)
    @test basis == MomentumBasis(3)
    @test hash(basis) == hash(MomentumBasis(3))

    momentum = @inferred LinearMomentum([1, -2, 0])
    @test momentum[1] == 1 // 1
    @test momentum[2] == -2 // 1
    @test !iszero(momentum)
    @test iszero(LinearMomentum([0, 0]))
    @test momentum == LinearMomentum([1, -2, 0])
    @test hash(momentum) == hash(LinearMomentum([1, -2, 0]))

    e1 = @inferred basis_momentum(basis, 1)
    e2 = @inferred basis_momentum(basis, 2)
    @test e1 + 2 * e2 == LinearMomentum([1, 2, 0])
    @test 3 * e1 - e2 == LinearMomentum([3, -1, 0])
    @test zero(momentum) == LinearMomentum([0, 0, 0])
    @test_throws DimensionMismatch LinearMomentum([1]) + LinearMomentum([1, 2])

    @test routing_recursively_concrete(typeof(basis))
    @test routing_recursively_concrete(typeof(momentum))
end

@testset "exact conservation-matrix routing" begin
    triangle = [
        -1 0 1
        1 -1 0
        0 1 -1
    ]
    routing = @inferred exact_momentum_routing(triangle)
    coefficients = @inferred routing_matrix(routing)

    @test routing isa MomentumRouting
    @test isempty(routing.edges)
    @test length(routing) == 3
    @test length(routing.basis) == 1
    @test routing.pivot_columns == [1, 2]
    @test routing.free_columns == [3]
    @test coefficients == reshape(Rational{Int}[1, 1, 1], 3, 1)
    @test triangle * coefficients == zeros(Rational{Int}, 3, 1)
    @test routing_recursively_concrete(typeof(routing))

    permuted_rows = triangle[[3, 1, 2], :]
    permuted_routing = @inferred exact_momentum_routing(permuted_rows)
    @test permuted_routing == routing
    @test hash(permuted_routing) == hash(routing)
end

@testset "tree and multiple loops" begin
    tree_edges = RoutingEdge[RoutingEdge(1, 2), RoutingEdge(2, 3)]
    tree = @inferred momentum_routing(tree_edges)
    @test length(tree.basis) == connected_cycle_rank(tree_edges) == 0
    @test size(routing_matrix(tree)) == (2, 0)
    @test independent_incidence(tree_edges) * routing_matrix(tree) ==
        zeros(Rational{Int}, 3, 0)

    theta_edges = RoutingEdge[
        RoutingEdge(1, 2, 1), RoutingEdge(1, 2, 2), RoutingEdge(1, 2, 3)
    ]
    theta = @inferred momentum_routing(theta_edges)
    @test length(theta.basis) == connected_cycle_rank(theta_edges) == 2
    @test routing_matrix(theta) == Rational{Int}[-1 -1; 1 0; 0 1]
    @test independent_incidence(theta_edges) * routing_matrix(theta) ==
        zeros(Rational{Int}, 2, 2)

    square_diagonal_edges = RoutingEdge[
        RoutingEdge(1, 2),
        RoutingEdge(2, 3),
        RoutingEdge(3, 4),
        RoutingEdge(4, 1),
        RoutingEdge(1, 3),
    ]
    square_diagonal = @inferred momentum_routing(square_diagonal_edges)
    @test length(square_diagonal.basis) == connected_cycle_rank(square_diagonal_edges) == 2
    @test length(square_diagonal) == 5
    @test independent_incidence(square_diagonal_edges) * routing_matrix(square_diagonal) ==
        zeros(Rational{Int}, 4, 2)
    @test routing_recursively_concrete(typeof(square_diagonal))
end

@testset "routing is independent of input edge order" begin
    edges = RoutingEdge[
        RoutingEdge(3, 1), RoutingEdge(1, 2), RoutingEdge(2, 3), RoutingEdge(1, 3)
    ]
    expected = @inferred momentum_routing(edges)
    permuted = @inferred momentum_routing(edges[[4, 2, 1, 3]])

    @test expected == permuted
    @test hash(expected) == hash(permuted)

    parallel = RoutingEdge[RoutingEdge(1, 2, 3), RoutingEdge(1, 2, 1), RoutingEdge(1, 2, 2)]
    parallel_permuted = parallel[[2, 3, 1]]
    p1 = @inferred momentum_routing(parallel)
    p2 = @inferred momentum_routing(parallel_permuted)
    @test p1 == p2
    @test hash(p1) == hash(p2)

    @test_throws ArgumentError momentum_routing(
        RoutingEdge[RoutingEdge(1, 2), RoutingEdge(1, 2)]
    )
end

@testset "self loops are exact cycle variables" begin
    edges = RoutingEdge[RoutingEdge(1, 1), RoutingEdge(1, 2)]
    routing = @inferred momentum_routing(edges)
    @test length(routing.basis) == connected_cycle_rank(edges) == 1
    @test routing_matrix(routing) == reshape(Rational{Int}[1, 0], 2, 1)
    @test independent_incidence(edges) * routing_matrix(routing) ==
        zeros(Rational{Int}, 2, 1)
end
