using KeldyshContraction, Test
using KeldyshContraction: Bulk, In, Out, Edge

@qfields ϕᶜ::Boson(Classical) ϕᴾ::Boson(Quantum)

@testset "construction" begin
    using KeldyshContraction: Diagram, Diagrams, Contraction, wick_contraction
    @inferred Diagrams{3,0}()

    @qfields c::Boson(Classical) q::Boson(Quantum)
    vs = KeldyshContraction.Contraction[
        (c(Out()), bar(q)),
        (c, bar(q)),
        (c, bar(q)(In())),
    ]
    @inferred Diagram(vs, Val(3), Val(0))
    @test_throws MethodError Diagram(vs)
end

@testset "static order entry point" begin
    using KeldyshContraction:
        In, Out, DressedPropagator, SelfEnergy, Diagram, FixedVector, matrix

    @qfields c::Boson(Classical) q::Boson(Quantum)
    L = InteractionLagrangian(
        -(1//2 * (c^2 + q^2) * bar(c) * bar(q) + 1//2 * c * q * (bar(c)^2 + bar(q)^2))
    )
    inout = c(Out()) * bar(q)(In())

    diagrams = @inferred wick_contraction(
        inout, L, Val(1), Val(3); simplify=false, _set_reg_to_zero=true
    )
    @test diagrams isa Diagrams{3,0}
    vacuum_diagrams = @inferred KeldyshContraction._wick_contraction(L.lagrangian, Val(2))
    @test vacuum_diagrams isa Diagrams{2,0}
    @test_throws MethodError KeldyshContraction._wick_contraction(L.lagrangian)

    propagator = @inferred DressedPropagator(
        L, Val(1), Val(3); simplify=false, _set_reg_to_zero=true
    )
    @test propagator isa DressedPropagator{3,0}
    @test @inferred(matrix(propagator)) isa Matrix{Diagrams{3,0}}
    self_energy = @inferred SelfEnergy(propagator, Val(1))
    @test self_energy isa SelfEnergy{1,0}
    @test @inferred(matrix(self_energy)) isa Matrix{Diagrams{1,0}}
    @test_throws AssertionError DressedPropagator(L, Val(1), Val(2))
    @test_throws AssertionError wick_contraction(inout, L, Val(1), Val(2))
    @test_throws AssertionError SelfEnergy(propagator, Val(2))
    @test @inferred(topologies(propagator.keldysh)) isa
        Dict{FixedVector{0,Int},Vector{Diagram{3,0}}}
    @test @inferred(wigner_transform(propagator)) isa DressedPropagator{3,0}

    @test_throws MethodError wick_contraction(inout, L, 1)
    @test_throws MethodError wick_contraction(inout, L, Val(1))
    @test_throws MethodError DressedPropagator(L, 1)
    @test_throws MethodError DressedPropagator(L, Val(1))
    @test_throws MethodError SelfEnergy(propagator, 1)
    @test_throws MethodError SelfEnergy(propagator)
end

@testset "prefactor multiplication" begin
    using KeldyshContraction: Diagram, Diagrams, Contraction

    @qfields c::Boson(Classical) q::Boson(Quantum)
    vs = KeldyshContraction.Contraction[
        (c(Out()), bar(q)),
        (c, bar(q)),
        (c, bar(q)(In())),
    ]
    d = Diagram(vs, Val(3), Val(0))
    ds = Diagrams(Dict(d => Complex{Rational{Int64}}(1.0)))
    ds2 = ds * 2.0
    ds2′ = 2.0 * ds
    @test isequal(ds2, ds2′)
end

@testset "is_connected" begin
    @qfields c::Boson(Classical) q::Boson(Quantum)

    vs = KeldyshContraction.Contraction[
        (c(Out()), bar(q)),
        (c, bar(q)),
        (c, bar(q)(In())),
    ]
    @test KeldyshContraction.is_connected(vs)

    vs2 = KeldyshContraction.Contraction[(c, bar(q))]
    @test KeldyshContraction.is_connected(vs2)

    vs3 = KeldyshContraction.Contraction[(c, bar(q)), (c(Out()), bar(q)(In()))]
    @test !KeldyshContraction.is_connected(vs3)
end

@testset "bulk multiplicity" begin
    using SmallCollections

    vs = FixedVector([(Out(), Bulk()), (Bulk(), Bulk()), (Bulk(), In())])
    vs = map(tt -> KeldyshContraction.index.(tt), vs)
    @test KeldyshContraction.bulk_multiplicity(vs) == Int[]

    vs2 = FixedVector([
        (Out(), Bulk()),
        (Bulk(), Bulk()),
        (Bulk(), Bulk(2)),
        (Bulk(2), Bulk(2)),
        (Bulk(2), In()),
    ])
    vs2 = map(tt -> KeldyshContraction.index.(tt), vs2)
    @test KeldyshContraction.bulk_multiplicity(vs2) == Int[1]

    vs3 = FixedVector([
        (Out(), Bulk()),
        (Bulk(), Bulk(2)),
        (Bulk(2), Bulk()),
        (Bulk(2), Bulk(2)),
        (Bulk(), In()),
    ])
    vs3 = map(tt -> KeldyshContraction.index.(tt), vs3)
    @test KeldyshContraction.bulk_multiplicity(vs3) == Int[2]

    vs4 = FixedVector([
        (Out(), Bulk()),
        (Bulk(1), Bulk(2)),
        (Bulk(2), Bulk(1)),
        (Bulk(2), Bulk(2)),
        (Bulk(1), Bulk(3)),
        (Bulk(3), Bulk(3)),
        (Bulk(3), In()),
    ])
    vs4 = map(tt -> KeldyshContraction.index.(tt), vs4)
    @test KeldyshContraction.bulk_multiplicity(vs4) == Int[2, 1, 0]
end

@testset "vertices" begin
    edges1 = Tuple{Int,Int}[]
    @test KeldyshContraction.vertices(edges1) == Set{Int}()

    edges2 = [(1, 2)]
    @test KeldyshContraction.vertices(edges2) == Set([1, 2])

    edges3 = [(1, 2), (2, 3), (3, 4)]
    @test KeldyshContraction.vertices(edges3) == Set([1, 2, 3, 4])

    edges4 = [(1, 1), (2, 2)]
    @test KeldyshContraction.vertices(edges4) == Set([1, 2])

    edges5 = [(1, 2), (2, 3), (1, 3)]
    @test KeldyshContraction.vertices(edges5) == Set([1, 2, 3])
end

@testset "connected components" begin
    edges1 = Tuple{Int,Int}[]
    vertices1 = Set{Int}()
    @test KeldyshContraction.connected_components(vertices1, edges1) == Vector{Set{Int}}()

    vertices2 = Set([1])
    edges2 = Tuple{Int,Int}[]
    comps2 = KeldyshContraction.connected_components(vertices2, edges2)
    @test length(comps2) == 1
    @test comps2[1] == Set([1])

    edges3 = [(1, 2), (2, 3), (3, 4)]
    vertices3 = KeldyshContraction.vertices(edges3)
    comps3 = KeldyshContraction.connected_components(vertices3, edges3)
    @test length(comps3) == 1
    @test comps3[1] == Set([1, 2, 3, 4])

    edges4 = [(1, 2), (3, 4), (5, 6)]
    vertices4 = KeldyshContraction.vertices(edges4)
    comps4 = KeldyshContraction.connected_components(vertices4, edges4)
    @test length(comps4) == 3
    @test Set([1, 2]) ∈ comps4
    @test Set([3, 4]) ∈ comps4
    @test Set([5, 6]) ∈ comps4

    edges5 = [(1, 2), (2, 3), (3, 1), (4, 5), (5, 6), (6, 4)]
    vertices5 = KeldyshContraction.vertices(edges5)
    comps5 = KeldyshContraction.connected_components(vertices5, edges5)
    @test length(comps5) == 2
    @test Set([1, 2, 3]) ∈ comps5
    @test Set([4, 5, 6]) ∈ comps5

    edges6 = [(1, 2), (2, 3)]
    vertices6 = Set([1, 2, 3, 4, 5])
    comps6 = KeldyshContraction.connected_components(vertices6, edges6)
    @test length(comps6) == 3
    @test Set([1, 2, 3]) ∈ comps6
    @test Set([4]) ∈ comps6
    @test Set([5]) ∈ comps6
end

@testset "Diagrams unique collection and prefactor sum" begin
    using KeldyshContraction: Diagram, Diagrams, Contraction
    c1 = (ϕᴾ, bar(ϕᶜ)(In()))
    c2 = (ϕᶜ, bar(ϕᶜ))
    c3 = (ϕᶜ(Out()), bar(ϕᴾ))
    contractions1 = Contraction[c1, c2, c3]
    contractions2 = Contraction[c1, c2, c3]
    contractions3 = Contraction[c2, c3, c1]
    contractions4 = Contraction[c1, c3]

    d1 = Diagram(contractions1, Val(3), Val(0))
    d2 = Diagram(contractions2, Val(3), Val(0))
    d3 = Diagram(contractions3, Val(3), Val(0))
    d4 = Diagram(contractions4, Val(2), Val(0))

    diagrams = Diagrams{3,0}()
    push!(diagrams, d1, 1.0)
    push!(diagrams, d2, 1.0)
    push!(diagrams, d3, 1.0)
    @test_throws MethodError push!(diagrams, d4, 1.0)

    collected = collect(diagrams)
    @test length(collected) == 1
    for (d, pref) in diagrams.diagrams
        if length(d.contractions) == 3
            @test pref == 3.0
        elseif length(d.contractions) == 2
            @test pref == 1.0
        end
    end
end
