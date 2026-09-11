using KeldyshContraction
using SymbolicUtils
using Test

const KC = KeldyshContraction

@testset "isequal/hash contract across symbolic representations" begin
    @qfields ϕ::Boson
    field = ϕ[Classical]

    scalar_one = one(field)
    scalar_zero = zero(field)
    single_field_mul = KC.QMul(1, Field{Boson}[field])
    single_field_add = KC.QAdd(KC.QMul{Int,Boson}[single_field_mul])

    for (a, b) in (
        (scalar_one, 1),
        (scalar_zero, 0),
        (single_field_mul, field),
        (single_field_add, field),
    )
        @test isequal(a, b)
        @test isequal(b, a)
        @test hash(a) == hash(b)
    end

    # The contract must also hold in hashed collections.
    @test Dict{Any,Int}(scalar_one => 1)[1] == 1
    @test Dict{Any,Int}(single_field_mul => 1)[field] == 1
    @test Set{Any}([single_field_add]) == Set{Any}([field])
end

@testset "QAdd owns its storage" begin
    @qfields ϕ::Boson ψ::Boson
    ϕc, ψq = ϕ[Classical], ψ[Quantum]

    caller_owned = KC.QMul{Int,Boson}[
        KC.QMul(1, Field{Boson}[ϕc]), KC.QMul(1, Field{Boson}[ψq])
    ]
    q = KC.QAdd(caller_owned)

    # Public construction must not retain the caller's mutable vector.
    empty!(caller_owned)
    @test length(KC.terms(q)) == 2

    # Package-native and SymbolicUtils accessors must not expose internal storage either.
    ts = KC.terms(q)
    empty!(ts)
    @test length(KC.terms(q)) == 2

    ast_args = SymbolicUtils.arguments(q)
    empty!(ast_args)
    @test length(KC.terms(q)) == 2
end
