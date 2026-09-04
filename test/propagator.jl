using KeldyshContraction, Test
using KeldyshContraction: In, Out, Classical, Quantum, is_in
using KeldyshContraction: Edge, position, contour
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus

@qfields c::Boson(Classical) q::Boson(Quantum)

@testset "propagator checks" begin
    @test_throws AssertionError Edge((c, c))
    @test_throws AssertionError Edge((c(In()), bar(c)(Out())))
    @test_throws AssertionError Edge((c(Out()), bar(c)(In())))
    @test_throws AssertionError Edge((c(Out()), bar(c)(Out())))
    @test_throws AssertionError Edge((q, bar(q)(In())))
    @test_throws AssertionError Edge((c(In()), bar(q)))
    @test_throws AssertionError Edge((c, bar(c)(Out())))
end

@testset "properties" begin
    using KeldyshContraction: PropagatorType, Regularisation
    p = Edge(q, bar(c)(In()))
    @test KeldyshContraction.position_category(p) == :in
    @test_broken KeldyshContraction.contours(p) == [Quantum, Classical]
    @test !KeldyshContraction.is_bulk(p)
    @test_broken KeldyshContraction.regularisations(p) == fill(Regularisation.Zero, 2)
    @test KeldyshContraction.propagator_type(p) == PropagatorType.Advanced
end

@testset "diagram construction" begin
    using KeldyshContraction: Diagram, Contraction
    contractions = Contraction[
        (q, bar(c)(In())),
        (c, bar(c)),
        (c(Out()), bar(q)),
    ]

    @inferred Diagram(contractions, Val(3), Val(0))
end

@testset "sort" begin
    using KeldyshContraction: Bulk, In, Out, sort_by_position_and_type
    p1 = (q, bar(c)(In()))
    p2 = (q, bar(c))
    @test isequal(sort!([p1, p2]; by=sort_by_position_and_type), [p2, p1])

    b1 = Bulk(1)
    b2 = Bulk(2)
    test1 = [
        (c(Out()), bar(q)(In())),
        (c(b1), bar(q)(b2)),
        (c(b1), bar(c)(b2)),
        (c(b2), bar(c)(b1)),
        (c(b2), bar(c)(In())),
    ]
    test2 = [
        (c(Out()), bar(q)(In())),
        (c(b1), bar(c)(b2)),
        (c(b1), bar(q)(b2)),
        (c(b2), bar(c)(b1)),
        (c(b2), bar(c)(In())),
    ]
    perm1 = sortperm(test1; by=sort_by_position_and_type)
    perm2 = sortperm(test2; by=sort_by_position_and_type)
    @test isequal(test1[perm1], test2[perm2])
    @test isequal(test2, test2[perm2])
end

@testset "adjoint" begin
    p1 = Edge(q, bar(c)(In()))
    p2 = Edge(c, bar(q)(In()))

    @test isequal(p1', p2)

    p = Edge(c, bar(c)(In()))
    @test_broken isequal(p', -1 * p)
end

@testset "regularisation" begin
    p = Edge(q(Plus), bar(c))
    @test_broken KeldyshContraction.regular(p) == false

    p = Edge(c(Minus), bar(c)(In()))
    @test_broken KeldyshContraction.regular(p) == true
end

@testset "propagator type" begin
    using KeldyshContraction: is_keldysh, is_retarded, is_advanced
    using KeldyshContraction: PropagatorType

    @test is_keldysh(PropagatorType.Keldysh)
    @test is_retarded(PropagatorType.Retarded)
    @test is_advanced(PropagatorType.Advanced)
    @test is_keldysh(Edge(c, bar(c)))
    @test is_retarded(Edge(c, bar(q)))
    @test is_advanced(Edge(q, bar(c)))
end

@testset "position" begin
    using KeldyshContraction: position, same_position, Bulk, In, Out
    p = (q, bar(c))
    @test same_position(p)
    p = (q(Out()), bar(c))
    @test !same_position(p)
    p = (q(Bulk(3)), bar(c)(Bulk(3)))
    @test same_position(p)
end

@testset "make spectral" begin
    using KeldyshContraction:
        make_spectral, make_retarded, make_advanced, PropagatorType, Bulk
    p = Edge(q, bar(c)(Bulk(2)))
    k = Edge(c, bar(q)(Bulk(2)))
    @test isequal(repr(make_spectral(p)), "A(y₁,y₂)")
    @test isequal(make_spectral(p), Edge(q, bar(c)(Bulk(2)), PropagatorType.Spectral))
    @test isequal(make_retarded(p), k)
    @test isequal(make_advanced(k), p)
end
