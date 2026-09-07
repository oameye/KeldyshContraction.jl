using KeldyshContraction, Test
using KeldyshContraction: matrix
import KeldyshContraction as KC

@qfields sum_ψ_raw::Fermion sum_χ_raw::Fermion
const sum_ψ = sum_ψ_raw
const sum_χ = sum_χ_raw
const sum_ψ₁ = sum_ψ[One]
const sum_ψ₂ = sum_ψ[Two]
const sum_χ₂ = sum_χ[Two]

function single_family_fermion_sum(coefficient::C) where {C<:Real}
    vertex = coefficient * sum_ψ₁ * sum_ψ₂ * bar(sum_ψ₁) * bar(sum_ψ₂)
    L_u = InteractionLagrangian(vertex, :u)
    L_v = InteractionLagrangian(vertex, :v)
    return L_u + L_v
end

function multi_family_fermion_sum(coefficient::C) where {C<:Real}
    vertex = coefficient * sum_ψ₁ * sum_χ₂ * bar(sum_ψ₁) * bar(sum_χ₂)
    L_u = InteractionLagrangian(vertex, :u)
    L_v = InteractionLagrangian(vertex, :v)
    return L_u + L_v
end

function same_components(a::DressedPropagator, b::DressedPropagator)
    return a.keldysh == b.keldysh && a.retarded == b.retarded && a.advanced == b.advanced
end

@testset "fermionic LagrangianSum exact propagation" begin
    Ls = @inferred single_family_fermion_sum(1 // 1)
    @test Ls isa KC.LagrangianSum{Rational{Int64},Fermion}

    G = @inferred DressedPropagator(Ls, Val(1), Val(3); simplify=false)
    u = parameter_monomial(:u)
    v = parameter_monomial(:v)
    @test Set(parameters(G)) == Set((u, v))
    @test G[u] isa DressedPropagator{KC.ComplexRationals,Fermion,1,3,0}
    @test same_components(G[u], G[v])

    Σ = @inferred SelfEnergy(G)
    @test Set(parameters(Σ)) == Set((u, v))
    @test Σ[u] isa SelfEnergy{KC.ComplexRationals,Fermion,1,1,0}
    @test iszero((@inferred matrix(G[u]))[2, 1])
    @test iszero((@inferred matrix(Σ[u]))[2, 1])
end

@testset "fermionic LagrangianSum mixed second order" begin
    Ls = @inferred single_family_fermion_sum(1 // 1)
    G = @inferred DressedPropagator(Ls, Val(2), Val(5); simplify=false)

    u = parameter_monomial(:u)
    v = parameter_monomial(:v)
    u², uv, v² = u^2, u * v, v^2
    @test Set(parameters(G)) == Set((u², uv, v²))

    G_u² = G[u²]
    G_uv = G[uv]
    G_v² = G[v²]
    @test same_components(G_u², G_v²)
    @test G_uv.keldysh == 2 * G_u².keldysh
    @test G_uv.retarded == 2 * G_u².retarded
    @test G_uv.advanced == 2 * G_u².advanced
    @test any(x -> !iszero(x), (G_u².keldysh, G_u².retarded, G_u².advanced))

    Σ = @inferred SelfEnergy(G)
    @test Set(parameters(Σ)) == Set((u², uv, v²))
end

@testset "fermionic LagrangianSum multi-family target" begin
    Ls = @inferred multi_family_fermion_sum(1 // 1)
    @test_throws ArgumentError DressedPropagator(Ls, Val(1), Val(3); simplify=false)

    G = @inferred DressedPropagator(Ls, Val(1), Val(3); target=sum_ψ, simplify=false)
    ps = parameters(G)
    @test Set(ps) == Set((parameter_monomial(:u), parameter_monomial(:v)))
    @test all(p -> G[p] isa DressedPropagator{KC.ComplexRationals,Fermion,1,3,0}, ps)
end

@testset "fermionic LagrangianSum floating coefficients" begin
    Ls = @inferred single_family_fermion_sum(1.0)
    G = @inferred DressedPropagator(Ls, Val(1), Val(3); simplify=false)
    Σ = @inferred SelfEnergy(G)

    @test G[:u] isa DressedPropagator{ComplexF64,Fermion,1,3,0}
    @test Σ[:u] isa SelfEnergy{ComplexF64,Fermion,1,1,0}
end
