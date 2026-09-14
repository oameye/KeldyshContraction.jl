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
end
