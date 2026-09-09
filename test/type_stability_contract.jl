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

@qfields contract_ψ_raw::Fermion contract_χ_raw::Fermion
const contract_ψ = contract_ψ_raw
const contract_χ = contract_χ_raw
const contract_ψ₁ = contract_ψ[One]
const contract_ψ₂ = contract_ψ[Two]
const contract_χ₂ = contract_χ[Two]

function assert_supported_public_contract(
    coefficient::C, ::Type{D}
) where {C<:Real,D<:Number}
    c = contract_c
    q = contract_q
    interaction =
        -coefficient * ((c^2 + q^2) * bar(c) * bar(q) + c * q * (bar(c)^2 + bar(q)^2))

    @test c isa Field{Boson}
    @test q isa Field{Boson}
    @test @inferred(field_family(c)) === contract_ϕ
    @test @inferred(field_family(q)) === contract_ϕ

    ∂xc = @inferred partial(c, :x)
    @test ∂xc isa Field{Boson}
    @test @inferred(field_family(∂xc)) === contract_ϕ
    @test @inferred(derivatives(∂xc)) == [:x]
    @test contract_recursively_concrete(typeof(∂xc))

    converted = @inferred convert_coefficients(D, interaction)
    rationalized = @inferred rationalize_coefficients(interaction)
    @test contract_recursively_concrete(typeof(converted))
    @test contract_recursively_concrete(typeof(rationalized))

    L = @inferred InteractionLagrangian(interaction)
    @test @inferred(field_families(L)) == [contract_ϕ]
    @test @inferred(target_family(L)) === contract_ϕ
    @test @inferred(parameters(L)) isa ParameterMonomial

    inout = c(Out()) * bar(q)(In())
    diagrams = @inferred(
        wick_contraction(inout, L, Val(1), Val(3); simplify=false, _set_reg_to_zero=true)
    )
    @test contract_recursively_concrete(typeof(diagrams))
    diagram_topologies = @inferred topologies(diagrams)
    @test contract_recursively_concrete(typeof(diagram_topologies))

    G = @inferred(DressedPropagator(L, Val(1), Val(3); simplify=false))
    Σ = @inferred SelfEnergy(G)
    Gk = @inferred fourier_transform(G)
    Σk = @inferred SelfEnergy(Gk)
    GW = @inferred wigner_transform(Gk; gradient_order=Val(0))
    ΣW = @inferred wigner_transform(Σk; gradient_order=Val(0))

    @test typeof(G).parameters[1] === D
    @test typeof(Σ).parameters[1] === D
    @test Gk isa FourierDressedPropagator{D,Boson,1,3,0}
    @test Σk isa FourierSelfEnergy{D,Boson,1,1,0}
    @test GW isa WignerDressedPropagator{D,Boson,1,3,0,0,HomogeneousWignerContext}
    @test ΣW isa WignerSelfEnergy{D,Boson,1,1,0,0,HomogeneousWignerContext}
    @test @inferred(parameters(G)) == parameters(L)
    @test @inferred(parameters(Σ)) == parameters(G)
    @test @inferred(parameters(Gk)) == parameters(G)
    @test @inferred(parameters(Σk)) == parameters(Gk)
    @test @inferred(parameters(GW)) == parameters(Gk)
    @test @inferred(parameters(ΣW)) == parameters(Σk)
    @test @inferred(gradient_order(GW)) == Val(0)
    @test @inferred(gradient_order(ΣW)) == Val(0)

    Gm = @inferred matrix(G)
    Σm = @inferred matrix(Σ)
    Gkm = @inferred matrix(Gk)
    Σkm = @inferred matrix(Σk)
    GWm = @inferred matrix(GW)
    ΣWm = @inferred matrix(ΣW)
    @test contract_recursively_concrete(typeof(Gm))
    @test contract_recursively_concrete(typeof(Σm))
    @test contract_recursively_concrete(typeof(Gkm))
    @test contract_recursively_concrete(typeof(Σkm))
    @test contract_recursively_concrete(typeof(GWm))
    @test contract_recursively_concrete(typeof(ΣWm))
    matrix_topologies = @inferred topologies(Gm[1, 1])
    @test contract_recursively_concrete(typeof(matrix_topologies))

    @test contract_recursively_concrete(typeof(L))
    @test contract_recursively_concrete(typeof(G))
    @test contract_recursively_concrete(typeof(Σ))
    @test contract_recursively_concrete(typeof(Gk))
    @test contract_recursively_concrete(typeof(Σk))
    @test contract_recursively_concrete(typeof(GW))
    @test contract_recursively_concrete(typeof(ΣW))

    io = IOBuffer()
    @test @inferred(show(io, c)) === nothing
    io = IOBuffer()
    @test @inferred(show(io, MIME"text/latex"(), c)) === nothing
    io = IOBuffer()
    @test @inferred(show(io, MIME"text/latex"(), L)) === nothing

    return nothing
end

function assert_supported_fermion_contract(
    coefficient::C, ::Type{D}
) where {C<:Real,D<:Number}
    ψ₁ = contract_ψ₁
    χ₂ = contract_χ₂
    interaction = coefficient * ψ₁ * χ₂ * bar(ψ₁) * bar(χ₂)

    @test ψ₁ isa Field{Fermion}
    @test χ₂ isa Field{Fermion}
    @test @inferred(field_family(ψ₁)) === contract_ψ
    @test @inferred(field_family(χ₂)) === contract_χ
    @test @inferred(bar(ψ₁)) isa Field{Fermion}

    ∂xχ₂ = @inferred partial(χ₂, :x)
    @test ∂xχ₂ isa Field{Fermion}
    @test @inferred(field_family(∂xχ₂)) === contract_χ
    @test @inferred(derivatives(∂xχ₂)) == [:x]
    @test contract_recursively_concrete(typeof(∂xχ₂))

    converted = @inferred convert_coefficients(D, interaction)
    rationalized = @inferred rationalize_coefficients(interaction)
    @test contract_recursively_concrete(typeof(converted))
    @test contract_recursively_concrete(typeof(rationalized))

    L = @inferred InteractionLagrangian(interaction)
    families = @inferred field_families(L)
    @test Set(families) == Set((contract_ψ, contract_χ))
    @test @inferred(target_family(L, contract_ψ)) === contract_ψ
    @test @inferred(parameters(L)) isa ParameterMonomial

    inout = ψ₁(Out()) * bar(ψ₁)(In())
    diagrams = @inferred(
        wick_contraction(inout, L, Val(1), Val(3); simplify=false, _set_reg_to_zero=true)
    )
    @test contract_recursively_concrete(typeof(diagrams))

    G = @inferred(DressedPropagator(L, Val(1), Val(3); target=contract_ψ, simplify=false))
    Σ = @inferred SelfEnergy(G)

    @test typeof(G).parameters[1] === D
    @test typeof(Σ).parameters[1] === D
    @test @inferred(parameters(G)) == parameters(L)
    @test @inferred(parameters(Σ)) == parameters(G)

    Gm = @inferred matrix(G)
    Σm = @inferred matrix(Σ)
    @test contract_recursively_concrete(typeof(Gm))
    @test contract_recursively_concrete(typeof(Σm))
    @test iszero(Gm[2, 1])
    @test iszero(Σm[2, 1])

    L_peer = @inferred InteractionLagrangian(interaction, :h)
    Ls = @inferred L + L_peer
    Gs = @inferred(DressedPropagator(Ls, Val(1), Val(3); target=contract_ψ, simplify=false))
    Σs = @inferred SelfEnergy(Gs)
    sum_parameters = @inferred parameters(Gs)
    self_energy_parameters = @inferred parameters(Σs)
    @test Set(sum_parameters) == Set((parameter_monomial(:g), parameter_monomial(:h)))
    @test Set(self_energy_parameters) == Set(sum_parameters)
    @test typeof(Gs[:g]).parameters[1] === D
    @test typeof(Σs[:g]).parameters[1] === D
    @test contract_recursively_concrete(typeof(Ls))
    @test contract_recursively_concrete(typeof(Gs))
    @test contract_recursively_concrete(typeof(Σs))

    derivative_interaction = coefficient * ψ₁ * ∂xχ₂ * bar(ψ₁) * bar(∂xχ₂)
    Ld = @inferred InteractionLagrangian(derivative_interaction, :d)
    Gd = @inferred(DressedPropagator(Ld, Val(1), Val(3); target=contract_ψ, simplify=false))
    Σd = @inferred SelfEnergy(Gd)
    Gdk = @inferred fourier_transform(Gd)
    Σdk = @inferred SelfEnergy(Gdk)
    GdW = @inferred wigner_transform(Gdk; gradient_order=Val(0))
    ΣdW = @inferred wigner_transform(Σdk; gradient_order=Val(0))
    @test typeof(Gd).parameters[1] === D
    @test typeof(Σd).parameters[1] === D
    @test Gdk isa FourierDressedPropagator{D,Fermion,1,3,0}
    @test Σdk isa FourierSelfEnergy{D,Fermion,1,1,0}
    @test GdW isa WignerDressedPropagator{D,Fermion,1,3,0,0,HomogeneousWignerContext}
    @test ΣdW isa WignerSelfEnergy{D,Fermion,1,1,0,0,HomogeneousWignerContext}
    @test contract_recursively_concrete(typeof(Ld))
    @test contract_recursively_concrete(typeof(Gd))
    @test contract_recursively_concrete(typeof(Σd))
    @test contract_recursively_concrete(typeof(Gdk))
    @test contract_recursively_concrete(typeof(Σdk))
    @test contract_recursively_concrete(typeof(GdW))
    @test contract_recursively_concrete(typeof(ΣdW))
    @test contract_recursively_concrete(typeof(@inferred matrix(GdW)))
    @test contract_recursively_concrete(typeof(@inferred matrix(ΣdW)))

    @test contract_recursively_concrete(typeof(L))
    @test contract_recursively_concrete(typeof(G))
    @test contract_recursively_concrete(typeof(Σ))

    return nothing
end

@testset "supported public type-stability contract" begin
    @testset "bosonic exact coefficients" begin
        @test @inferred(assert_supported_public_contract(1 // 2, KC.ComplexRationals)) ===
            nothing
    end

    @testset "bosonic floating-point coefficients" begin
        @test @inferred(assert_supported_public_contract(0.5, ComplexF64)) === nothing
    end

    @testset "fermionic exact coefficients" begin
        @test @inferred(assert_supported_fermion_contract(1 // 2, KC.ComplexRationals)) ===
            nothing
    end

    @testset "fermionic floating-point coefficients" begin
        @test @inferred(assert_supported_fermion_contract(0.5, ComplexF64)) === nothing
    end
end
