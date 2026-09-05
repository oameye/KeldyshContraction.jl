using KeldyshContraction, Test
using KeldyshContraction: In, Out, Bulk
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus
import KeldyshContraction as KC

@testset "concrete Field{Boson}" begin
    c = @inferred Field{Boson}(:c, Classical)
    q = @inferred Field{Boson}(:q, Quantum)

    @test typeof(c) === Field{Boson}
    @test isconcretetype(typeof(c))
    @test KC.keldysh_index(c) === Classical
    @test KC.keldysh_index(q) === Quantum
    @test KC.is_classical(c)
    @test KC.is_quantum(q)
    @test KC.is_unbarred(c)

    cbar = @inferred bar(c)
    @test typeof(cbar) === Field{Boson}
    @test KC.is_barred(cbar)
    @test KC.keldysh_index(cbar) === Classical
    @test KC.position(cbar) == KC.position(c)
    @test KC.regularisation(cbar) == KC.regularisation(c)
    @test @inferred(bar(cbar)) == c

    @qfields ϕ::Boson(Classical) ψ::Boson(Quantum)
    @test typeof(ϕ) === Field{Boson}
    @test typeof(ψ) === Field{Boson}
end

@testset "concrete field metadata" begin
    using KeldyshContraction: FieldIndex, FieldIndices, field_indices

    spin = FieldIndex(:spin, 1)
    flavor = FieldIndex(:flavor, 2)
    indices = @inferred FieldIndices(spin, flavor)
    f = @inferred Field{Boson}(
        :ψ, Classical, KC.Orientation.Unbarred, KC.Regularisation.Zero, Bulk(), indices
    )

    @test isconcretetype(typeof(indices))
    @test field_indices(f) == indices
    @test collect(field_indices(f)) == [spin, flavor]
    @test hash(f) == hash(f)
end

@testset "concrete QMul" begin
    @qfields c::Boson(Classical) q::Boson(Quantum)

    mul = @inferred KC.QMul(1, Field{Boson}[c, q])
    @test typeof(mul) === KC.QMul{Int,Boson}
    @test eltype(mul.args_nc) === Field{Boson}
    @test isconcretetype(typeof(mul))

    @test typeof(@inferred(c * q)) === KC.QMul{Int,Boson}
    @test typeof(@inferred(0.5 * c * q)) === KC.QMul{Float64,Boson}
    @test typeof(@inferred(c^2)) === KC.QMul{Int,Boson}

    # Bosonic canonical ordering is value-independent and sign-free.
    @test isequal(c * q, q * c)
    @test isequal(bar(c) * c, c * bar(c))
end

@testset "closed symbolic zero and one" begin
    @qfields c::Boson(Classical) q::Boson(Quantum)

    z = @inferred zero(c)
    o = @inferred one(c)
    @test typeof(z) === KC.QMul{Int,Boson}
    @test typeof(o) === KC.QMul{Int,Boson}
    @test iszero(z)
    @test isone(o)
    @test isempty(z.args_nc)
    @test isempty(o.args_nc)

    @test typeof(@inferred(c^0)) === KC.QMul{Int,Boson}
    @test isone(c^0)
    @test typeof(@inferred(c + 0)) === KC.QAdd{Int,Boson}
    @test isequal(c + 0, KC.QAdd(Field{Boson}[c]))

    zero_mul = @inferred KC.QMul(0.0, Field{Boson}[c, q])
    @test typeof(zero_mul) === KC.QMul{Float64,Boson}
    @test iszero(zero_mul)
    @test isempty(zero_mul.args_nc)
end

@testset "concrete QAdd and coefficient promotion" begin
    using KeldyshContraction: QAdd, QMul
    @qfields c::Boson(Classical) q::Boson(Quantum)

    a = @inferred QAdd(Field{Boson}[c, q])
    @test typeof(a) === QAdd{Int,Boson}
    @test eltype(a.arguments) === QMul{Int,Boson}
    @test isconcretetype(typeof(a))

    mixed = @inferred(2.0 * c * c + 2 * q * q)
    @test typeof(mixed) === QAdd{Float64,Boson}
    @test eltype(mixed.arguments) === QMul{Float64,Boson}

    complex_mixed = @inferred(0.5 * c * c + 2im * q * q)
    @test typeof(complex_mixed) === QAdd{ComplexF64,Boson}

    @test typeof(@inferred(a * a)) === QAdd{Int,Boson}
    @test typeof(@inferred(a + c)) === QAdd{Int,Boson}
end

@testset "explicit coefficient conversion" begin
    using KeldyshContraction: QAdd, QMul
    @qfields c::Boson(Classical) q::Boson(Quantum)

    expr = 0.5 * (c^2 + q^2) * bar(c) * bar(q)
    @test typeof(expr) === QAdd{Float64,Boson}

    rational_expr = @inferred rationalize_coefficients(expr)
    @test typeof(rational_expr) === QAdd{Rational{Int64},Boson}
    @test typeof(expr) === QAdd{Float64,Boson}

    converted = @inferred convert_coefficients(Float32, expr)
    @test typeof(converted) === QAdd{Float32,Boson}

    L = InteractionLagrangian(expr)
    @test typeof(L.lagrangian) === QAdd{Float64,Boson}
end

@testset "SymbolicUtils / TermInterface interop" begin
    using TermInterface, SymbolicUtils
    @qfields ϕ::Boson(Classical) ψ::Boson(Quantum)

    @test TermInterface.head(ϕ) == :call
    @test !SymbolicUtils.iscall(ϕ)
    @test SymbolicUtils.iscall(ϕ * ψ)
    @test SymbolicUtils.iscall(ϕ + ψ)
    @test SymbolicUtils.operation(ϕ + ψ) == +
    @test SymbolicUtils.operation(ϕ * ψ) == *

    mul = 2 * ϕ * ψ
    args = @inferred SymbolicUtils.arguments(mul)
    @test eltype(args) === Union{Int,Field{Boson}}
    @test first(args) == 2
    @test args[2:end] == KC.fields(mul)
    @test isnothing(TermInterface.metadata(mul))
    @test isequal(@inferred(TermInterface.maketerm(typeof(mul), *, args, nothing)), mul)

    unit_mul = ϕ * ψ
    @test isequal(
        @inferred(TermInterface.maketerm(typeof(unit_mul), *, Field{Boson}[ϕ, ψ], nothing)),
        unit_mul,
    )
end

@testset "position and regularisation are value data" begin
    @qfields ϕ::Boson(Classical) ψ::Boson(Quantum)

    @test typeof(@inferred(ϕ(Bulk(3)))) === Field{Boson}
    @test typeof(@inferred(ψ(In()))) === Field{Boson}
    @test typeof(@inferred(ψ(Out()))) === Field{Boson}
    @test typeof(@inferred(ϕ(Plus))) === Field{Boson}
    @test typeof(@inferred(ϕ(Minus))) === Field{Boson}

    @test KC.position(ψ(In())) == In()
    @test KC.position(ψ(Out())) == Out()
    @test KC.regularisation(ϕ(Plus)) === KC.Regularisation.Plus

    to_sort = [ϕ(Bulk(3)), ϕ(Bulk(1)), ϕ, ψ(In()), ψ(Out())]
    @test eltype(to_sort) === Field{Boson}
    @test @inferred KC.is_bulk(ϕ)
    @test @inferred KC.is_bulk(ϕ * ψ)
end

@testset "bosonic interaction expressions" begin
    @qfields ϕ::Boson(Classical) ψ::Boson(Quantum)

    elastic = -0.5 * (bar(ϕ) * bar(ψ) * (ϕ^2 + ψ^2) + (bar(ϕ)^2 + bar(ψ)^2) * ϕ * ψ)
    loss =
        0.5 * bar(ϕ) * bar(ψ) * (ϕ(Minus)^2 + ψ(Minus)^2) -
        0.5 * ϕ(Plus) * ψ(Plus) * (bar(ϕ)^2 + bar(ψ)^2) +
        bar(ϕ) * bar(ψ) * (ϕ(Plus)^2 + ϕ(Minus)^2)

    @test elastic isa KC.QAdd{Float64,Boson}
    @test loss isa KC.QAdd{Float64,Boson}
    @test @inferred KC.is_conserved(first(elastic.arguments))
    @test @inferred KC.is_physical(first(elastic.arguments))
end
