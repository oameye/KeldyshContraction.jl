using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields ψ::Fermion

@testset "fermionic Keldysh labels" begin
    ψ₁ = ψ[One]
    ψ₂ = ψ[Two]

    @test ψ isa FieldFamily{Fermion}
    @test ψ₁ isa Field{Fermion}
    @test ψ₂ isa Field{Fermion}
    @test One === KC.KeldyshIndex.First
    @test Two === KC.KeldyshIndex.Second
    @test KC.keldysh_index(ψ₁) === One
    @test KC.keldysh_index(ψ₂) === Two
    @test KC.is_one(ψ₁)
    @test KC.is_two(ψ₂)
    @test !KC.is_one(ψ₂)
    @test !KC.is_two(ψ₁)

    barred = bar(ψ₁)
    @test barred isa Field{Fermion}
    @test field_family(barred) === ψ
    @test KC.keldysh_index(barred) === One
    @test KC.is_barred(barred)
    @test bar(barred) == ψ₁

    @test @inferred(field_family(ψ₁)) === ψ
    @test @inferred(bar(ψ₁)) isa Field{Fermion}
end

@testset "fermionic Grassmann algebra" begin
    ψ₁ = ψ[One]
    ψ₂ = ψ[Two]

    forward = @inferred ψ₁ * ψ₂
    reverse = @inferred ψ₂ * ψ₁
    @test forward == -reverse
    @test typeof(forward) === typeof(reverse)

    nilpotent = @inferred ψ₁ * ψ₁
    @test iszero(nilpotent)
    @test typeof(nilpotent) === typeof(forward)
    @test isempty(KC.fields(nilpotent))

    squared = @inferred ψ₁^2
    @test iszero(squared)
    @test typeof(squared) === typeof(forward)

    barred_partner = @inferred ψ₁ * bar(ψ₁)
    @test !iszero(barred_partner)
    @test typeof(barred_partner) === typeof(forward)

    different_component = @inferred ψ₁ * ψ₂
    @test !iszero(different_component)

    left_associated = @inferred((ψ₁ * ψ₂ * bar(ψ₁)) * bar(ψ₂))
    product_reference = @inferred((ψ₁ * ψ₂) * (bar(ψ₁) * bar(ψ₂)))
    @test left_associated == product_reference

    sum = ψ₁ + ψ₂
    right_multiplied = @inferred sum * bar(ψ₁)
    distributed = @inferred ψ₁ * bar(ψ₁) + ψ₂ * bar(ψ₁)
    @test right_multiplied == distributed
end

function reference_permutation_sign(permutation)
    inversions = 0
    for i in 1:(length(permutation) - 1), j in (i + 1):length(permutation)
        inversions += permutation[i] > permutation[j]
    end
    return iseven(inversions) ? Int8(1) : Int8(-1)
end

@testset "fermionic Wick permutation parity" begin
    for permutation in ([1], [1, 2], [2, 1], [1, 3, 2], [2, 3, 1], [3, 2, 1])
        @test KC.pairing_sign(Fermion, permutation) ==
            reference_permutation_sign(permutation)
    end
end

@testset "fermionic LO propagator mapping" begin
    ψ₁ = ψ[One]
    ψ₂ = ψ[Two]
    bψ₁ = bar(ψ₁)
    bψ₂ = bar(ψ₂)

    @test KC.propagator_type(ψ₁, bψ₁) === KC.PropagatorType.Retarded
    @test KC.propagator_type(ψ₁, bψ₂) === KC.PropagatorType.Keldysh
    @test KC.propagator_type(ψ₂, bψ₂) === KC.PropagatorType.Advanced
    @test_throws ArgumentError KC.propagator_type(ψ₂, bψ₁)

    @test KC.contraction_filter(KC.Contraction(ψ₁, bψ₁))
    @test KC.contraction_filter(KC.Contraction(ψ₁, bψ₂))
    @test KC.contraction_filter(KC.Contraction(ψ₂, bψ₂))
    @test !KC.contraction_filter(KC.Contraction(ψ₂, bψ₁))
    @test_throws AssertionError KC.Edge(ψ₂, bψ₁)

    advanced = KC.Edge(ψ₂, bψ₂)
    retarded = adjoint(advanced)
    @test KC.is_advanced(advanced)
    @test KC.is_retarded(retarded)
    @test KC.is_one(retarded.out)
    @test KC.is_one(retarded.in)
end

@testset "fermionic LO matrices" begin
    ψ₁ = ψ[One]
    ψ₂ = ψ[Two]

    dR = KC.Diagram([KC.Edge(ψ₁, bar(ψ₁))], Val(1), Val(0))
    dK = KC.Diagram([KC.Edge(ψ₁, bar(ψ₂))], Val(1), Val(0))
    dA = KC.Diagram([KC.Edge(ψ₂, bar(ψ₂))], Val(1), Val(0))

    R = KC.Diagrams([dR], ComplexF64(1))
    K = KC.Diagrams([dK], ComplexF64(2))
    A = KC.Diagrams([dA], ComplexF64(3))
    parameter = parameter_monomial(:g)

    G = DressedPropagator(K, R, A, Val(1), parameter)
    Gm = @inferred KC.matrix(G)
    @test Gm[1, 1] == R
    @test Gm[1, 2] == K
    @test iszero(Gm[2, 1])
    @test Gm[2, 2] == A
    @test eltype(Gm) === KC.Diagrams{ComplexF64,Fermion,1,0}

    Σ = SelfEnergy{ComplexF64,Fermion,1,1,0}(K, R, A, parameter)
    Σm = @inferred KC.matrix(Σ)
    @test Σm[1, 1] == R
    @test Σm[1, 2] == K
    @test iszero(Σm[2, 1])
    @test Σm[2, 2] == A
    @test eltype(Σm) === KC.Diagrams{ComplexF64,Fermion,1,0}
end

@testset "fermionic multi-family first-order reference" begin
    @qfields ψref::Fermion χref::Fermion
    ψ₁ = ψref[One]
    χ₂ = χref[Two]

    vertex = @inferred ψ₁ * χ₂ * bar(ψ₁) * bar(χ₂)
    L = @inferred InteractionLagrangian(vertex, :u)

    @test Set(field_families(L)) == Set((ψref, χref))
    @test_throws ArgumentError DressedPropagator(L, Val(1), Val(3); simplify=false)

    G = @inferred DressedPropagator(L, Val(1), Val(3); target=ψref, simplify=false)
    @test G isa DressedPropagator{KC.ComplexRationals,Fermion,1,3,0}
    @test length(G.retarded) == 1
    @test length(G.keldysh) == 1
    @test isempty(G.advanced)
    @test only(values(G.retarded.diagrams)) == -im
    @test only(values(G.keldysh.diagrams)) == -im

    retarded_edges = KC.contractions(only(keys(G.retarded.diagrams)))
    @test count(KC.is_retarded, retarded_edges) == 2
    @test count(KC.is_advanced, retarded_edges) == 1
    @test count(KC.is_keldysh, retarded_edges) == 0

    keldysh_edges = KC.contractions(only(keys(G.keldysh.diagrams)))
    @test count(KC.is_retarded, keldysh_edges) == 1
    @test count(KC.is_advanced, keldysh_edges) == 1
    @test count(KC.is_keldysh, keldysh_edges) == 1

    Σ = @inferred SelfEnergy(G)
    @test Σ isa SelfEnergy{KC.ComplexRationals,Fermion,1,1,0}
    @test length(Σ.retarded) == 1
    @test isempty(Σ.keldysh)
    @test isempty(Σ.advanced)
    @test only(values(Σ.retarded.diagrams)) == -im

    internal_edge = only(KC.contractions(only(keys(Σ.retarded.diagrams))))
    @test KC.is_advanced(internal_edge)
    @test field_family(internal_edge.out) === χref
    @test field_family(internal_edge.in) === χref

    Gm = @inferred KC.matrix(G)
    Σm = @inferred KC.matrix(Σ)
    @test iszero(Gm[2, 1])
    @test iszero(Σm[2, 1])
end
