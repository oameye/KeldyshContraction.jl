using KeldyshContraction, Test
using KeldyshContraction: In, Out, Bulk
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus
import KeldyshContraction as KC

@testset "typestable QSym" begin
    @inferred Destroy Destroy(:ψ, Classical)
    @inferred Create Create(:ψ, Classical)

    @qfields c::Destroy(Classical) q::Destroy(Quantum)
end

@testset "Type stable QMul" begin
    @qfields c::Destroy(Classical) q::Destroy(Quantum)

    @inferred KC.QMul(1, [c, c])
    @inferred c * c
    @inferred KC.QMul(1, [c, c]) * KC.QMul(1.0, [c, c])
    @inferred c^2
    @inferred 1.0 * c

    mul = c * c
    @inferred c * mul
    @inferred mul * mul

    @inferred 0.5 * q^2 * c' * q'
    example = 0.5 * q * c * c' * q'
    # @inferred arguments(example)
    # example = 0.5 * q*c * c' * q'
    # @code_warntype InteractionLagrangian(example)
    # ^ Does not have to be type stable, as it is called only once
end

@testset "Type stable QAdd" begin
    @qfields c::Destroy(Classical) q::Destroy(Quantum)

    @inferred KC.QAdd([c, c])
    @inferred c + c
    @inferred -c
    @inferred 2.0 * c * c
    @inferred 2.0 * c * c + 2 * c * c
    @inferred c * c + 0

    add = c + c
    @inferred add + c
    @inferred add + 0.0

    @inferred 0.5 * (c^2 + q^2) * c' * q'
    @inferred 0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2)
    @inferred 0.5 * (c^2 + q^2) * c' * q' + 2 * c * q * ((c')^2 + (q')^2)
    @inferred 0.5 * (c^2 + q^2) * c' * q' + 2 * im * c * q * ((c')^2 + (q')^2)
end

@testset "SymbolicUtils interface" begin
    using TermInterface, SymbolicUtils
    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)

    @test TermInterface.head(ϕ) == :call
    @test SymbolicUtils.iscall(ϕ) == false
    @test SymbolicUtils.iscall(ϕ * ϕ) == true
    @test SymbolicUtils.iscall(ϕ + ϕ) == true
    @test SymbolicUtils.operation(ϕ + ϕ) == +
    @test SymbolicUtils.operation(ϕ * ϕ) == *
    @test SymbolicUtils.arguments(2 * ϕ * ϕ) == [2, ϕ, ϕ]
    @test isnothing(TermInterface.metadata(2 * ϕ * ϕ))
    @test isequal(TermInterface.maketerm(KC.QMul, *, [ϕ, ϕ], nothing), ϕ * ϕ)

    @testset "SymbolicUtils promotion" begin
        # Test the promotion of Keldysh fields
        typeof(ϕ + 0.0)
        @test SymbolicUtils.promote_symtype(+, ϕ, 0.0) isa KC.QField broken = true
        @test SymbolicUtils.promote_symtype(*, ϕ, 1) isa KC.QField broken = true
    end
end

@testset "normal ordering" begin
    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)
    @test isequal(ϕ' * ϕ, ϕ * ϕ')
end

@testset "ones and zeros" begin
    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)
    @test isone(ϕ) == false
    @test iszero(ϕ) == false
    @test one(ϕ) == 1
    @test zero(ϕ) == 0
end

# @testset "hash" begin # how to test hashes?
#     @test hash(ϕ * ϕ, hash(4)) == 0x85ebad55106e875a
#     @test hash(ϕ + ϕ, hash(4)) == 0x8fc9006919ba7ef2
# end

@testset "isequal" begin
    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)
    # Test the equality of two Keldysh fields
    @test ϕ == ϕ
    @test isequal(ϕ * ϕ, ϕ * ϕ)
    @test isequal(ϕ * ψ, ϕ * ψ)
    @test isequal(ϕ + ψ, ϕ + ψ)

    @test isequal(ψ + ϕ, ϕ + ψ) broken = true
    @test isequal(ψ * ϕ, ϕ * ψ) broken = true
    @test isequal(0.0 + ϕ, ϕ + 0)
    ϕ2 = ϕ + ϕ
    @test isequal(ϕ2 + 0, ϕ + ϕ + 0)
    @test isequal(ϕ2 + ϕ, ϕ + ϕ + ϕ)
end

@testset "simplification" begin
    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)
    using SymbolicUtils
    @test isequal(ϕ + ϕ, 2 * ϕ) broken = true
    @test isequal(ϕ + ϕ + ϕ, 3 * ϕ) broken = true

    @test isequal((ϕ + ϕ) * (ϕ + ϕ), 4 * ϕ^2) broken = true
    @test isequal((ϕ + ϕ) * (ϕ + ϕ), ϕ^2 + ϕ^2 + ϕ^2 + ϕ^2)
    # SymbolicUtils.simplify((ϕ + ϕ) * (ψ + ϕ) + 3 * (ϕ + ϕ) * (ψ + ϕ))
    # SymbolicUtils.expand((ϕ + ϕ) * (ψ + ϕ) + 3 * (ϕ + ϕ) * (ψ + ϕ))
    # ^ broken due to giving args::Vector{Any}
    # check if it is still happens after QAdd being type stable
end

@testset "adjoint" begin
    using KeldyshContraction: is_creation, is_annihilation, is_conserved

    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)
    ϕ′ = Create(KC.name(ϕ), Classical, KC.regularisation(ϕ), KC.position(ϕ))
    ψ′ = Create(KC.name(ψ), Quantum, KC.regularisation(ψ), KC.position(ψ))

    @test is_creation(ϕ′)
    @test !is_creation(ϕ)
    @test is_annihilation(ϕ)
    @test !is_conserved(ϕ)

    # Test the adjoint of Keldysh fields
    @test isequal(ϕ', ϕ′)
    @test adjoint(ψ) == ψ′
    @test adjoint(ψ′) == ψ
    @test isequal(adjoint(ϕ * ψ), ϕ′ * ψ′)
    @test isequal(adjoint(ϕ + ψ), ϕ′ + ψ′)
end

@testset "position" begin
    using KeldyshContraction: position, is_bulk
    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)
    # Test the position function
    @test KC.position(ϕ).index == 1
    @test KC.position(ψ(In())) == In()
    @test KC.position(ψ(Out())) == Out()
    # @test KC.position(ϕ + ψ) == [0]
    # @test KC.position(ϕ * ψ) == [0]

    # @test KC.position(ϕ + ψ(In)) == [0, 1]
    # @test KC.position(ϕ * ψ(In)) == [0, 1]
    to_sort = [ϕ(Bulk(3)), ϕ(Bulk(1)), ϕ, ψ(In()), ψ(Out())]
    sorted = [ψ(Out()), ϕ, ϕ(Bulk(1)), ϕ(Bulk(3)), ψ(In())]
    @test isequal(sort(to_sort; by=KC.position), sorted)

    @inferred is_bulk(ϕ)
    @inferred is_bulk(ϕ * ψ)
end

@testset "more math" begin
    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)

    # Test the math operations
    @test isequal(ϕ / 2, 0.5 * ϕ)
    @test isequal(ϕ//2, 0.5 * ϕ) skip = true

    @test isequal((ϕ^2), ϕ * ϕ)

    @test isequal(ϕ, ϕ + 0)
    @test isequal(0 + ϕ, ϕ)

    # @code_warntype isequal(ϕ, ϕ + 0)

    mul = ϕ * ϕ
    add = ϕ + ϕ
    @test isequal(ϕ * mul, ϕ^3)
    @test isequal(mul * ϕ, ϕ^3)
    @test isequal(ϕ + add, ϕ + ϕ + ϕ)
    @test isequal(add + ϕ, ϕ + ϕ + ϕ)
    @test isequal(mul * add, ϕ^3 + ϕ^3)
    @test isequal(add * mul, ϕ^3 + ϕ^3)
    @test isequal(add + mul, ϕ + ϕ + ϕ^2)
    @test isequal(mul + add, ϕ + ϕ + ϕ^2)
end

@testset "quantum-classical" begin
    using KeldyshContraction: is_quantum, is_classical
    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)

    @test is_quantum(ψ)
    @test is_classical(ϕ)
    @test !is_quantum(ϕ)
    @test !is_classical(ψ)
end

@testset "some lagrangians" begin
    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)
    eleactic2boson = -0.5 * (ϕ' * ψ' * (ϕ^2 + ψ^2) + (ϕ'^2 + ψ'^2) * ϕ * ψ)

    loss2boson =
        0.5 * ϕ' * ψ' * (ϕ(Minus)^2 + ψ(Minus)^2) -
        0.5 * ϕ(Plus) * ψ(Plus) * (ϕ'^2 + ψ'^2) + ϕ' * ψ' * (ϕ(Plus)^2 + ϕ(Minus)^2)
end

@testset "is_conserved" begin
    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)
    using KeldyshContraction: is_conserved

    @test !is_conserved(KeldyshContraction.QSym[])
    @test !is_conserved(KeldyshContraction.QSym[ϕ])

    @inferred is_conserved(ϕ * ψ)
    # @code_warntype is_conserved(ϕ * ψ)
end

@testset "is_physical" begin
    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)
    using KeldyshContraction: is_physical

    @inferred is_physical([ϕ, ψ])
    @inferred is_physical(ϕ * ψ)
end

@testset "QMul" begin
    using KeldyshContraction: QMul
    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)

    @test isequal(ϕ, QMul(1, [ϕ]))
    @test isequal(QMul(1, [ϕ]), ϕ)

    @test !isequal(ϕ, QMul(1, [ϕ, ϕ]))
    @test !isequal(QMul(1, [ϕ, ϕ]), ϕ)

    @test isequal(QMul(0, [ϕ]), 0)
    @test isequal(0, QMul(0, [ϕ]))

    @test !isequal(QMul(0, [ϕ]), 1)
    @test !isequal(1, QMul(0, [ϕ]))

    @test iszero(QMul(0, [ϕ]))
    @test iszero(zero(ϕ * ϕ))

    # type promotion
    promote_type(KC.QMul{Int64}, KC.QMul{Float64})
end

@testset "conversion/promotion" begin
    using KeldyshContraction: QAdd, QMul

    @qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)

    a = QAdd([ϕ, ψ])
    b = QAdd([QMul(1, [ϕ]), QMul(1, [ψ])])
    @test isequal(a, b)
    @test typeof(a + b) == QAdd{Float64}
    @test typeof(a * b) == QAdd{Float64}
    @test typeof(convert(QAdd{Float64}, a)) == QAdd{Float64}
    @test typeof(convert(QMul{Float64}, QMul(1, [ϕ]))) == QMul{Float64}
end
