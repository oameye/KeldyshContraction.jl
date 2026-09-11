using KeldyshContraction, Test
import KeldyshContraction as KC

function recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        recursively_concrete(eltype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        recursively_concrete(FT, seen) || return false
    end
    return true
end

@testset "static computational spine" begin
    @qfields ϕ::Boson
    c, q = ϕ[Classical], ϕ[Quantum]
    elastic = -0.5 * ((c^2 + q^2) * bar(c) * bar(q) + c * q * (bar(c)^2 + bar(q)^2))

    L = @inferred InteractionLagrangian(elastic)
    G = @inferred DressedPropagator(L, Val(1), Val(3); simplify=false)
    Σ = @inferred SelfEnergy(G)
    Gk = @inferred fourier_transform(G)
    Σk = @inferred SelfEnergy(Gk)
    collision = @inferred KC.CollisionIntegral(Σ)

    GT = typeof(G)
    ΣT = typeof(Σ)
    @test GT.parameters[2] === Boson
    @test GT.parameters[3] == 1
    @test GT.parameters[4] == 3
    @test GT.parameters[5] == 0
    @test ΣT.parameters[2] === Boson
    @test ΣT.parameters[3] == 1
    @test ΣT.parameters[4] == 1
    @test ΣT.parameters[5] == 0

    @test Gk isa FourierDressedPropagator{ComplexF64,Boson,1,3,0}
    @test Σk isa FourierSelfEnergy{ComplexF64,Boson,1,1,0}
    @test recursively_concrete(typeof(L))
    @test recursively_concrete(typeof(G))
    @test recursively_concrete(typeof(Σ))
    @test recursively_concrete(typeof(Gk))
    @test recursively_concrete(typeof(Σk))
    @test recursively_concrete(typeof(collision))
    @test parameters(G) isa ParameterMonomial
    @test parameters(Σ) == parameters(G)
    @test parameters(Gk) == parameters(G)
    @test parameters(Σk) == parameters(Gk)

    Gm = @inferred KC.matrix(G)
    @test Gm[1, 1] === G.keldysh
    @test Gm[1, 2] === G.retarded
    @test Gm[2, 1] === G.advanced
    @test iszero(Gm[2, 2])

    Σm = @inferred KC.matrix(Σ)
    @test iszero(Σm[1, 1])
    @test Σm[1, 2] === Σ.advanced
    @test Σm[2, 1] === Σ.retarded
    @test Σm[2, 2] === Σ.keldysh

    Gkm = @inferred KC.matrix(Gk)
    @test Gkm[1, 1] === Gk.keldysh
    @test Gkm[1, 2] === Gk.retarded
    @test Gkm[2, 1] === Gk.advanced
    @test iszero(Gkm[2, 2])

    Σkm = @inferred KC.matrix(Σk)
    @test iszero(Σkm[1, 1])
    @test Σkm[1, 2] === Σk.advanced
    @test Σkm[2, 1] === Σk.retarded
    @test Σkm[2, 2] === Σk.keldysh
end
