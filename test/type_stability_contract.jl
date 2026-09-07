using KeldyshContraction, Test
using KeldyshContraction: In, Out, matrix
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

function assert_supported_public_contract(coefficient::C, ::Type{D}) where {C<:Real,D<:Number}
    c = contract_c
    q = contract_q
    interaction =
        -coefficient * ((c^2 + q^2) * bar(c) * bar(q) + c * q * (bar(c)^2 + bar(q)^2))

    @test c isa Field{Boson}
    @test q isa Field{Boson}
    @test @inferred(field_family(c)) === contract_ϕ
    @test @inferred(field_family(q)) === contract_ϕ

    converted = @inferred convert_coefficients(D, interaction)
    rationalized = @inferred rationalize_coefficients(interaction)
    @test contract_recursively_concrete(typeof(converted))
    @test contract_recursively_concrete(typeof(rationalized))

    L = @inferred InteractionLagrangian(interaction)
    @test @inferred(field_families(L)) == [contract_ϕ]
    @test @inferred(target_family(L)) === contract_ϕ
    @test @inferred(parameters(L)) isa ParameterMonomial

    inout = c(Out()) * bar(q)(In())
    diagrams = @inferred wick_contraction(
        inout, L, Val(1), Val(3); simplify=false, _set_reg_to_zero=true
    )
    @test contract_recursively_concrete(typeof(diagrams))
    @test contract_recursively_concrete(typeof(@inferred(topologies(diagrams)))

    G = @inferred DressedPropagator(L, Val(1), Val(3); simplify=false)
    Σ = @inferred SelfEnergy(G)
    Gk = @inferred wigner_transform(G)
    Σk = @inferred wigner_transform(Σ)

    @test typeof(G).parameters[1] === D
    @test typeof(Σ).parameters[1] === D
    @test typeof(Gk) === typeof(G)
    @test typeof(Σk) === typeof(Σ)
    @test @inferred(parameters(G)) == parameters(L)
    @test @inferred(parameters(Σ)) == parameters(G)

    Gm = @inferred matrix(G)
    Σm = @inferred matrix(Σ)
    @test contract_recursively_concrete(typeof(Gm))
    @test contract_recursively_concrete(typeof(Σm))
    @test contract_recursively_concrete(typeof(@inferred(topologies(Gm[1, 1])))

    @test contract_recursively_concrete(typeof(L))
    @test contract_recursively_concrete(typeof(G))
    @test contract_recursively_concrete(typeof(Σ))
    @test contract_recursively_concrete(typeof(Gk))
    @test contract_recursively_concrete(typeof(Σk))

    io = IOBuffer()
    @test @inferred(show(io, c)) === nothing
    io = IOBuffer()
    @test @inferred(show(io, MIME"text/latex"(), c)) === nothing
    io = IOBuffer()
    @test @inferred(show(io, MIME"text/latex"(), L)) === nothing

    return nothing
end

@testset "supported public type-stability contract" begin
    @testset "exact coefficients" begin
        @test @inferred(assert_supported_public_contract(1 // 2, KC.ComplexRationals)) ===
            nothing
    end

    @testset "floating-point coefficients" begin
        @test @inferred(assert_supported_public_contract(0.5, ComplexF64)) === nothing
    end
end
