using KeldyshContraction, Test
import KeldyshContraction as KC

function contract_recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        contract_recursively_concrete(eltype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        contract_recursively_concrete(FT, seen) || return false
    end
    return true
end

@qfields contract_ϕ::Boson
const contract_c = contract_ϕ[Classical]
const contract_q = contract_ϕ[Quantum]

function assert_supported_type_contract(coefficient::C, ::Type{D}) where {C<:Real,D<:Number}
    c = contract_c
    q = contract_q
    interaction =
        -coefficient * ((c^2 + q^2) * bar(c) * bar(q) + c * q * (bar(c)^2 + bar(q)^2))

    L = @inferred InteractionLagrangian(interaction)
    G = @inferred DressedPropagator(L, Val(1), Val(3); simplify=false)
    Σ = @inferred SelfEnergy(G)
    Gk = @inferred wigner_transform(G)
    Σk = @inferred wigner_transform(Σ)
    collision = @inferred KC.CollisionIntegral(Σk)

    @test typeof(G).parameters[1] === D
    @test typeof(Σ).parameters[1] === D
    @test typeof(Gk) === typeof(G)
    @test typeof(Σk) === typeof(Σ)
    @test collision isa KC.CollisionIntegral{D,0}

    @test contract_recursively_concrete(typeof(L))
    @test contract_recursively_concrete(typeof(G))
    @test contract_recursively_concrete(typeof(Σ))
    @test contract_recursively_concrete(typeof(Gk))
    @test contract_recursively_concrete(typeof(Σk))
    @test contract_recursively_concrete(typeof(collision))

    diagram = first(keys(Gk.keldysh.diagrams))
    io = IOBuffer()
    @test @inferred(show(io, diagram)) === nothing
    io = IOBuffer()
    @test @inferred(show(io, MIME"text/latex"(), diagram)) === nothing

    return nothing
end

@testset "supported type-stability contract" begin
    @testset "exact coefficients" begin
        @test @inferred(assert_supported_type_contract(1 // 2, KC.ComplexRationals)) ===
            nothing
    end

    @testset "floating-point coefficients" begin
        @test @inferred(assert_supported_type_contract(0.5, ComplexF64)) === nothing
    end
end
