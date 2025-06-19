using KeldyshContraction, Test
using KeldyshContraction: In, Out, Classical, Quantum, is_in
using KeldyshContraction: Edge, position, contour
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus

@qfields c::Destroy(Classical) q::Destroy(Quantum)

@testset "propagator checks" begin
    @test_throws MethodError Edge((c, c)) # annihilation creation
    @test_throws AssertionError Edge((c(In()), c'(Out()))) # In-Out
    @test_throws AssertionError Edge((c(Out()), c'(In()))) # In-Out
    @test_throws AssertionError Edge((c(Out()), c'(Out()))) # same coordinate
    @test_throws AssertionError Edge((q, q(In())')) # quantum-quantum
    @test_throws AssertionError Edge((c(In()), q')) # Out is In
    @test_throws AssertionError Edge((c, c'(Out()))) # In is Out
end

@testset "properties" begin
    using KeldyshContraction: PropagatorType, Regularisation
    p = Edge(q, c'(In()))
    @test KeldyshContraction.position_category(p) == :in
    @test_broken KeldyshContraction.contours(p) == [Quantum, Classical] # not needed
    @test !KeldyshContraction.is_bulk(p)
    @test_broken KeldyshContraction.regularisations(p) == fill(Regularisation.Zero, 2) # not needed
    @test KeldyshContraction.propagator_type(p) == PropagatorType.Advanced
end

@testset "diagram construction" begin
    using KeldyshContraction: Diagram, Contraction
    contractions = Contraction[(q, c'(In())), (c, c'), (c(Out()), q')]

    Diagram(contractions)
end

@testset "sort" begin
    using KeldyshContraction: Bulk, In, Out, sort_by_position_and_type
    p1 = (q, c'(In()))
    p2 = (q, c')
    @test isequal(sort!([p1, p2]; by=sort_by_position_and_type), [p2, p1])

    b1 = Bulk(1)
    b2 = Bulk(2)
    test1 = [
        (c(Out()), q'(In())),
        (c(b1), q'(b2)),
        (c(b1), c'(b2)),
        (c(b2), c'(b1)),
        (c(b2), c'(In())),
    ]
    test2 = [
        (c(Out()), q'(In())),
        (c(b1), c'(b2)),
        (c(b1), q'(b2)),
        (c(b2), c'(b1)),
        (c(b2), c'(In())),
    ]
    perm1 = sortperm(test1; by=sort_by_position_and_type)
    perm2 = sortperm(test2; by=sort_by_position_and_type)
    @test isequal(test1[perm1], test2[perm2])
    @test isequal(test2, test2[perm2])
end

@testset "adjoint" begin
    p1 = Edge(q, c'(In()))
    p2 = Edge(c, q'(In()))

    @test isequal(p1', p2) # (G^R)† = G^A

    p = Edge(c, c'(In()))
    @test_broken isequal(p', -1 * p) # (G^K)† = -G^K
end

@testset "regularisation" begin
    p = Edge(q(Plus), c')
    @test_broken KeldyshContraction.regular(p) == false

    p = Edge(c(Minus), c'(In()))
    @test_broken KeldyshContraction.regular(p) == true
end

@testset "propagator type" begin
    using KeldyshContraction: is_keldysh, is_retarded, is_advanced
    using KeldyshContraction: PropagatorType

    @test is_keldysh(PropagatorType.Keldysh)
    @test is_retarded(PropagatorType.Retarded)
    @test is_advanced(PropagatorType.Advanced)
    @test is_keldysh(Edge(c, c'))
    @test is_retarded(Edge(c, q'))
    @test is_advanced(Edge(q, c'))
end

@testset "position" begin
    using KeldyshContraction: position, same_position, Bulk, In, Out
    p = (q, c')
    @test same_position(p)
    p = (q(Out()), c')
    @test !same_position(p)
    p = (q(Bulk(3)), c(Bulk(3))')
    @test same_position(p)
end
