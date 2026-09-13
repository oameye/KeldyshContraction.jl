using KeldyshContraction, Test
using KeldyshContraction: In, Out, Bulk
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus
import KeldyshContraction as KC

@testset "concrete Field{Boson}" begin
    @qfields family::Boson
    c = @inferred family[Classical]
    q = @inferred family[Quantum]

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

    @test family isa FieldFamily{Boson}
    @test field_family(c) === family
    @test field_family(q) === family
end

@testset "concrete field metadata" begin
    using KeldyshContraction: FieldIndex, FieldIndices, field_indices

    spin = FieldIndex(:spin, 1)
    flavor = FieldIndex(:flavor, 2)
    indices = @inferred FieldIndices(spin, flavor)
    family = FieldFamily{Boson}(:ψ, indices)
    f = @inferred Field(
        family, Classical, KC.Orientation.Unbarred, KC.Regularisation.Zero, Bulk()
    )

    @test isconcretetype(typeof(indices))
    @test field_family(f) === family
    @test field_indices(f) == indices
    @test collect(field_indices(f)) == [spin, flavor]
    @test hash(f) == hash(f)
end

@testset "concrete QMul" begin
    @qfields family::Boson
    c, q = family[Classical], family[Quantum]

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
    @qfields family::Boson
    c, q = family[Classical], family[Quantum]

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
    @qfields family::Boson
    c, q = family[Classical], family[Quantum]

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
    @qfields family::Boson
    c, q = family[Classical], family[Quantum]

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
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

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
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

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
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

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

#####################################################################
# Recursive concreteness
#
# `isconcretetype` only inspects the outer type: it says nothing about what a
# `Vector` field actually stores. These walk the whole storage graph instead.
#####################################################################

"""Whether `T` and every type reachable through its fields are concrete."""
function is_recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        is_recursively_concrete(eltype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        is_recursively_concrete(FT, seen) || return false
    end
    return true
end

@testset "recursively concrete storage" begin
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

    @test is_recursively_concrete(Field{Boson})
    @test is_recursively_concrete(KC.FieldIndices)
    @test is_recursively_concrete(KC.Position)
    @test is_recursively_concrete(typeof(ϕ * ψ))
    @test is_recursively_concrete(typeof(ϕ + ψ))
    @test is_recursively_concrete(typeof(0.5 * (ϕ^2 + ψ^2) * bar(ϕ) * bar(ψ)))

    # The walker must be able to fail, or the assertions above prove nothing.
    @test !is_recursively_concrete(KC.QMul)
    @test !is_recursively_concrete(Vector{Any})
    @test !is_recursively_concrete(Tuple{Vector{KC.QField},Int})

    # Storage stays concrete no matter which coefficient type arithmetic lands on.
    for expr in (ϕ * ψ, 2 * ϕ, 0.5 * ϕ, 2im * ϕ, ϕ + ψ, (ϕ + ψ) * (ϕ + ψ), ϕ // 2)
        @test is_recursively_concrete(typeof(expr))
    end
end

@testset "statistics-dispatched exchange sign" begin
    using KeldyshContraction: exchange_sign, is_exchange_sign_free, canonicalize_fields!
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

    @test @inferred(exchange_sign(Boson)) === Int8(1)
    @test @inferred(exchange_sign(Boson, ϕ, ψ)) === Int8(1)
    @test @inferred(is_exchange_sign_free(Boson))

    # Canonicalization reports the sign it accumulated, and bosons never flip.
    unsorted = Field{Boson}[ψ, ϕ]
    @test @inferred(canonicalize_fields!(unsorted)) === Int8(1)
    @test issorted(unsorted)
    @test canonicalize_fields!(Field{Boson}[]) === Int8(1)

    # Bosonic ordering is sign-free, so reordering a product cannot change it.
    @test isequal(ϕ * ψ, ψ * ϕ)
    @test KC.coefficient(ϕ * ψ) == KC.coefficient(ψ * ϕ)
end

@testset "algebraic identities" begin
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

    @test isequal(ϕ / 2, 0.5 * ϕ)
    @test isequal(ϕ // 2, (1 // 2) * ϕ)
    @test isequal((ϕ * ψ) / 2, 0.5 * ϕ * ψ)
    @test isequal((ϕ + ψ) / 2, 0.5 * ϕ + 0.5 * ψ)
    @test isequal((ϕ + ψ) // 2, (1 // 2) * ϕ + (1 // 2) * ψ)

    @test isequal(ϕ^2, ϕ * ϕ)
    @test isequal((ϕ * ψ)^2, ϕ * ψ * ϕ * ψ)
    @test isequal((ϕ + ψ)^2, (ϕ + ψ) * (ϕ + ψ))
    @test isone(ϕ^0)
    # `n` must not be a literal: Julia rewrites `x^-1` to `inv(x)` before dispatch.
    n = -1
    @test_throws DomainError ϕ^n
    @test_throws DomainError (ϕ * ψ)^n
    @test_throws DomainError (ϕ + ψ)^n

    @test isequal(-ϕ, -1 * ϕ)
    @test isequal(-(ϕ * ψ), -1 * ϕ * ψ)
    @test isequal(-(ϕ + ψ), -1 * ϕ + -1 * ψ)
    @test isequal(ϕ - ψ, ϕ + (-ψ))
    @test isequal(2 - ϕ, 2 + (-ϕ))
    @test isequal(ϕ - 2, ϕ + (-2))

    @test isequal(ϕ, ϕ + 0)
    @test isequal(0 + ϕ, ϕ)

    mul, add = ϕ * ϕ, ϕ + ϕ
    @test isequal(ϕ * mul, ϕ^3)
    @test isequal(mul * ϕ, ϕ^3)
    @test isequal(ϕ + add, ϕ + ϕ + ϕ)
    @test isequal(add + ϕ, ϕ + ϕ + ϕ)
    @test isequal(mul * add, ϕ^3 + ϕ^3)
    @test isequal(add * mul, ϕ^3 + ϕ^3)
    @test isequal(add + mul, ϕ + ϕ + ϕ^2)
    @test isequal(mul + add, ϕ + ϕ + ϕ^2)
end

@testset "ones, zeros, equality and hashing" begin
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

    @test !isone(ϕ)
    @test !iszero(ϕ)
    # `one`/`zero` stay inside the IR rather than collapsing to a raw scalar.
    @test one(ϕ) isa KC.QMul{Int,Boson}
    @test zero(ϕ) isa KC.QMul{Int,Boson}
    @test one(ϕ) == 1
    @test zero(ϕ) == 0
    @test isone(one(ϕ))
    @test iszero(zero(ϕ))
    @test isone(one(Field{Boson}))
    @test iszero(zero(Field{Boson}))

    @test ϕ == ϕ
    @test isequal(ϕ * ϕ, ϕ * ϕ)
    @test isequal(ϕ * ψ, ϕ * ψ)
    @test isequal(ϕ + ψ, ϕ + ψ)
    @test isequal(ψ * ϕ, ϕ * ψ)
    @test isequal(0.0 + ϕ, ϕ + 0)

    # `QAdd` stores its summands in insertion order, so addition is not commutative.
    # Sorting them canonically would reorder the terms of every Lagrangian.
    @test isequal(ψ + ϕ, ϕ + ψ) broken = true

    ϕ2 = ϕ + ϕ
    @test isequal(ϕ2 + 0, ϕ + ϕ + 0)
    @test isequal(ϕ2 + ϕ, ϕ + ϕ + ϕ)

    @test hash(ϕ + ψ) == hash(ϕ + ψ)
    @test hash(ϕ * ψ) == hash(ϕ * ψ)
    @test hash(ϕ) == hash(ϕ)
    @test hash(bar(ϕ)) != hash(ϕ)
end

@testset "simplification leaves like terms uncollected" begin
    using SymbolicUtils
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

    # `QAdd` keeps every summand separately: like terms are never collected, and a sum
    # of one term never collapses to a `QMul`. Hence ϕ + ϕ stays a two-term `QAdd`.
    @test isequal(ϕ + ϕ, 2 * ϕ) broken = true
    @test isequal(ϕ + ϕ + ϕ, 3 * ϕ) broken = true
    @test isequal((ϕ + ϕ) * (ϕ + ϕ), 4 * ϕ^2) broken = true
    @test isequal((ϕ + ϕ) * (ϕ + ϕ), ϕ^2 + ϕ^2 + ϕ^2 + ϕ^2)

    @test length(SymbolicUtils.arguments(SymbolicUtils.expand((ϕ + ϕ) * (ψ + ϕ)))) == 4
    @test isequal(
        SymbolicUtils.simplify((ϕ + ϕ) * (ψ + ϕ) + 3 * (ϕ + ϕ) * (ψ + ϕ)),
        SymbolicUtils.expand((ϕ + ϕ) * (ψ + ϕ) + 3 * (ϕ + ϕ) * (ψ + ϕ)),
    )
end

@testset "bar is the orientation operation" begin
    using KeldyshContraction: is_quantum, is_classical, is_barred, is_unbarred
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

    @test is_quantum(ψ)
    @test is_classical(ϕ)
    @test !is_quantum(ϕ)
    @test !is_classical(ψ)

    # `bar` toggles orientation and preserves every other property.
    ϕ′ = bar(ϕ)
    @test is_barred(ϕ′) && is_unbarred(ϕ)
    @test KC.name(ϕ′) === KC.name(ϕ)
    @test KC.keldysh_index(ϕ′) === KC.keldysh_index(ϕ)
    @test KC.position(ϕ′) == KC.position(ϕ)
    @test KC.regularisation(ϕ′) === KC.regularisation(ϕ)
    @test KC.field_indices(ϕ′) == KC.field_indices(ϕ)
    @test field_family(ϕ′) === field_family(ϕ)
    @test isequal(bar(bar(ϕ)), ϕ)

    @test isequal(bar(ϕ * ψ), bar(ϕ) * bar(ψ))
    @test isequal(bar(ϕ + ψ), bar(ϕ) + bar(ψ))

    # `adjoint` is reserved for objects with a genuine mathematical adjoint.
    @test_throws MethodError ϕ'
end

@testset "conservation and physicality" begin
    using KeldyshContraction: is_conserved, is_physical
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

    @test !is_conserved(Field{Boson}[])
    @test !is_conserved(Field{Boson}[ϕ])
    @test !is_conserved(())
    @test !is_conserved(ϕ)
    @test is_conserved(Field{Boson}[ϕ, bar(ϕ)])
    @test is_conserved((ϕ, bar(ϕ)))
    @test @inferred is_conserved(ϕ * bar(ϕ))
    @test @inferred is_physical(Field{Boson}[ϕ, ψ])
    @test @inferred is_physical(ϕ * ψ)

    # An In() field must be barred and an Out() field unbarred.
    @test !is_physical(ϕ(In()))
    @test is_physical(bar(ϕ)(In()))
    @test is_physical(ϕ(Out()))
    @test !is_physical(bar(ϕ)(Out()))
end

@testset "QMul equality and promotion" begin
    using KeldyshContraction: QMul, QAdd
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

    @test isequal(ϕ, QMul(1, Field{Boson}[ϕ]))
    @test isequal(QMul(1, Field{Boson}[ϕ]), ϕ)
    @test !isequal(ϕ, QMul(1, Field{Boson}[ϕ, ϕ]))
    @test !isequal(QMul(1, Field{Boson}[ϕ, ϕ]), ϕ)

    @test isequal(QMul(0, Field{Boson}[ϕ]), 0)
    @test isequal(0, QMul(0, Field{Boson}[ϕ]))
    @test !isequal(QMul(0, Field{Boson}[ϕ]), 1)
    @test iszero(QMul(0, Field{Boson}[ϕ]))
    @test iszero(zero(ϕ * ϕ))

    @test promote_type(QMul{Int,Boson}, QMul{Float64,Boson}) === QMul{Float64,Boson}
    @test promote_type(QAdd{Int,Boson}, QAdd{Float64,Boson}) === QAdd{Float64,Boson}

    # The public vector constructor is defensive: it must not adopt the caller's vector.
    caller_owned = Field{Boson}[ψ, ϕ]
    mul = QMul(1, caller_owned)
    @test caller_owned == Field{Boson}[ψ, ϕ]
    @test issorted(KC.fields(mul))

    # Accessors hand back copies, so mutating them cannot corrupt the expression.
    fs = KC.fields(mul)
    empty!(fs)
    @test length(KC.fields(mul)) == 2
end

@testset "conversion and promotion" begin
    using KeldyshContraction: QAdd, QMul
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

    a = QAdd(Field{Boson}[ϕ, ψ])
    b = QAdd(QMul{Int,Boson}[QMul(1, Field{Boson}[ϕ]), QMul(1, Field{Boson}[ψ])])
    @test isequal(a, b)
    @test typeof(a + b) === QAdd{Int,Boson}
    @test typeof(a * b) === QAdd{Int,Boson}
    @test typeof(convert(QAdd{Float64,Boson}, a)) === QAdd{Float64,Boson}
    @test typeof(convert(QMul{Float64,Boson}, QMul(1, Field{Boson}[ϕ]))) ===
        QMul{Float64,Boson}
    @test convert(QAdd{Int,Boson}, a) === a

    # Coefficient promotion is decided by the argument types, never by their values.
    @test typeof(1 * ϕ + 1.0 * ψ) === QAdd{Float64,Boson}
    @test typeof(1.0 * ϕ + 1 * ψ) === QAdd{Float64,Boson}
    @test typeof(0.0 * ϕ + 1 * ψ) === QAdd{Float64,Boson}
    @test typeof((1 // 2) * ϕ * ψ) === QMul{Rational{Int},Boson}
end

@testset "explicit rationalisation" begin
    using KeldyshContraction: QAdd
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

    for coeff in (0.5, 0.3, 0.1)
        expr = coeff * (ϕ^2 + ψ^2) * bar(ϕ) * bar(ψ)
        @test typeof(rationalize_coefficients(expr)) === QAdd{Rational{Int},Boson}
    end

    mixed = 0.1 * (ϕ^2 + ψ^2) * bar(ϕ) * bar(ψ) + 0.5 * (ϕ^2 + ψ^2) * bar(ϕ) * bar(ψ)
    @test typeof(rationalize_coefficients(mixed)) === QAdd{Rational{Int},Boson}

    # Rationalising an already-exact expression is the identity on its type.
    exact = (1 // 2) * ϕ * ψ
    @test typeof(rationalize_coefficients(exact)) === typeof(exact)
end

@testset "SymbolicUtils promotion" begin
    using SymbolicUtils
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

    @test SymbolicUtils.promote_symtype(+, typeof(ϕ), Float64) <: KC.QField
    @test SymbolicUtils.promote_symtype(*, typeof(ϕ), Int) <: KC.QField
    @test SymbolicUtils.promote_symtype(*, Int, typeof(ϕ)) <: KC.QField
    @test SymbolicUtils.symtype(ϕ) === Field{Boson}
end

@testset "the IR has no empty sum" begin
    using KeldyshContraction: QAdd, QMul
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

    # An empty `QAdd` would make `iszero` true, `isone` false and `first(terms(q))` throw.
    @test_throws ArgumentError QAdd(QMul{Int,Boson}[])
    @test_throws ArgumentError QAdd{Int,Boson}(QMul{Int,Boson}[])

    # Algebraic zero is the one-term zero monomial, and stays inside the IR.
    z = zero(QAdd{Int,Boson})
    @test length(KC.terms(z)) == 1
    @test iszero(z)
    @test !isone(z)
    @test isone(one(QAdd{Int,Boson}))

    # Cancelling every term still yields a well-formed sum, not an empty one.
    cancelled = 0 * (ϕ + ψ)
    @test cancelled isa QAdd{Int,Boson}
    @test iszero(cancelled)
    @test length(KC.terms(cancelled)) == 1
end

@testset "field metadata rejects bad input" begin
    using KeldyshContraction: FieldIndex, FieldIndices, NO_FIELD_INDEX, MAX_FIELD_INDICES

    @test_throws ArgumentError FieldIndex(:charge, 1)
    @test_throws ArgumentError FieldIndices(FieldIndex(:spin, 1), FieldIndex(:spin, 2))
    @test_throws ArgumentError FieldIndices(FieldIndex(:band, NO_FIELD_INDEX))

    idx = FieldIndices(FieldIndex(:spin, 1), FieldIndex(:band, -3))
    @test length(idx) == 2
    @test eltype(idx) === FieldIndex
    @test collect(idx) isa Vector{FieldIndex}
    @test_throws BoundsError idx[0]
    @test_throws BoundsError idx[3]

    # Every supported kind round-trips, and ordering is canonical, not insertion order.
    full = FieldIndices(
        FieldIndex(:species, 4),
        FieldIndex(:band, 3),
        FieldIndex(:flavor, 2),
        FieldIndex(:spin, 1),
    )
    @test length(full) == MAX_FIELD_INDICES
    @test [i.value for i in full] == [1, 2, 3, 4]

    @test FieldIndices() < idx
    @test hash(idx) == hash(FieldIndices(FieldIndex(:band, -3), FieldIndex(:spin, 1)))
end

@testset "positions reject bad input" begin
    using KeldyshContraction: Position, subtraction, Regularisation

    @test_throws ArgumentError Bulk(0)
    @test_throws ArgumentError Bulk(-1)
    @test_throws ArgumentError subtraction(Regularisation.T[Regularisation.Plus])
    @test subtraction((Regularisation.Plus, Regularisation.Minus)) == 2
    @test KC.is_bulk(Bulk(2)) && !KC.is_bulk(In()) && !KC.is_bulk(Out())
    @test KC.swap_in_out(In()) == Out()
    @test KC.swap_in_out(Out()) == In()
    @test KC.swap_in_out(Bulk(2)) == Bulk(2)

    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]
    to_sort = [ϕ(Bulk(3)), ϕ(Bulk(1)), ϕ, ψ(In()), ψ(Out())]
    sorted = [ψ(Out()), ϕ, ϕ(Bulk(1)), ϕ(Bulk(3)), ψ(In())]
    @test isequal(sort(to_sort; by=KC.position), sorted)
end

@testset "printing covers every expression shape" begin
    using KeldyshContraction: QAdd, QMul, Momenta, Momentum
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

    # `QAdd` prints through the generic `QTerm` method; `QMul` has its own.
    @test occursin("+", repr(ϕ + ψ))
    @test occursin("*", repr(ϕ * ψ))
    @test !isempty(repr(bar(ϕ)))
    @test occursin("̄", repr(bar(ϕ)))

    # Regularisation superscripts.
    @test occursin("⁺", repr(ϕ(Plus)))
    @test occursin("⁻", repr(ϕ(Minus)))

    # A coefficient-only monomial has no fields left to print.
    @test repr(zero(ϕ)) == "0"
    @test repr(one(ϕ)) == "1"

    # No momenta at all prints as nothing; a stored zero prefactor prints as a bare "0".
    @test repr(Momenta()) == ""
    @test repr(Momenta([0], [Momentum(1)], Val(:raw))) == "0"
    @test repr(Momenta(1)) == repr(Momenta([1], [Momentum(1)]))

    @test !isempty(repr(2 * ϕ * ψ + 3 * ψ))
    @test !isempty(sprint(show, ϕ + ψ))
end

@testset "interop and scalar predicates on QField" begin
    using TermInterface, SymbolicUtils
    using KeldyshContraction: QField, QSym, QTerm, QMul, QAdd, Position
    @qfields family::Boson
    ϕ, ψ = family[Classical], family[Quantum]

    # Type-level `iscall` is used by SymbolicUtils rewriting.
    @test SymbolicUtils.iscall(QMul{Int,Boson})
    @test SymbolicUtils.iscall(QAdd{Int,Boson})
    @test isnothing(TermInterface.metadata(ϕ))

    # A bare field is neither the algebraic zero nor the algebraic one.
    @test !isone(ϕ)
    @test !iszero(ϕ)

    @test Int(Position(3)) == 3
    @test Int(KC.position(ϕ)) == 1
    @test Position(1) < Position(2)
end
