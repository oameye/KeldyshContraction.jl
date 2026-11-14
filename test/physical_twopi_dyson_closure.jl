using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields physical_dyson_ψ::Boson physical_dyson_χ::Boson

const PhysicalDysonCoeff = Complex{Rational{Int}}
const PhysicalDysonWord = Tuple{Vararg{Symbol}}

"Exact test-only polynomial in ordered noncommuting component symbols."
struct PhysicalDysonExpr
    terms::Dict{PhysicalDysonWord,PhysicalDysonCoeff}
end

PhysicalDysonExpr() = PhysicalDysonExpr(Dict{PhysicalDysonWord,PhysicalDysonCoeff}())
Base.zero(::Type{PhysicalDysonExpr}) = PhysicalDysonExpr()
Base.zero(::PhysicalDysonExpr) = PhysicalDysonExpr()
Base.iszero(expression::PhysicalDysonExpr) = isempty(expression.terms)
Base.isequal(a::PhysicalDysonExpr, b::PhysicalDysonExpr) = isequal(a.terms, b.terms)
Base.:(==)(a::PhysicalDysonExpr, b::PhysicalDysonExpr) = isequal(a, b)

function physical_dyson_scalar(value::Number)
    coefficient = convert(PhysicalDysonCoeff, value)
    terms = Dict{PhysicalDysonWord,PhysicalDysonCoeff}()
    iszero(coefficient) || (terms[()] = coefficient)
    return PhysicalDysonExpr(terms)
end

function physical_dyson_symbol(name::Symbol)
    return PhysicalDysonExpr(
        Dict{PhysicalDysonWord,PhysicalDysonCoeff}((name,) => one(PhysicalDysonCoeff))
    )
end

physical_dyson_expr(value::PhysicalDysonExpr) = value
physical_dyson_expr(value::Number) = physical_dyson_scalar(value)
physical_dyson_expr(value::Symbol) = physical_dyson_symbol(value)

function Base.:+(a::PhysicalDysonExpr, b::PhysicalDysonExpr)
    terms = copy(a.terms)
    for (word, coefficient) in b.terms
        combined = get(terms, word, zero(PhysicalDysonCoeff)) + coefficient
        if iszero(combined)
            delete!(terms, word)
        else
            terms[word] = combined
        end
    end
    return PhysicalDysonExpr(terms)
end

Base.:-(a::PhysicalDysonExpr) = (-1) * a
Base.:-(a::PhysicalDysonExpr, b::PhysicalDysonExpr) = a + (-b)

function Base.:*(factor::Number, expression::PhysicalDysonExpr)
    coefficient = convert(PhysicalDysonCoeff, factor)
    iszero(coefficient) && return PhysicalDysonExpr()
    return PhysicalDysonExpr(
        Dict{PhysicalDysonWord,PhysicalDysonCoeff}(
            word => coefficient * value for (word, value) in expression.terms
        ),
    )
end
Base.:*(expression::PhysicalDysonExpr, factor::Number) = factor * expression

function Base.:*(a::PhysicalDysonExpr, b::PhysicalDysonExpr)
    terms = Dict{PhysicalDysonWord,PhysicalDysonCoeff}()
    for (left_word, left_coefficient) in a.terms
        for (right_word, right_coefficient) in b.terms
            word = (left_word..., right_word...)
            combined =
                get(terms, word, zero(PhysicalDysonCoeff)) +
                left_coefficient * right_coefficient
            if iszero(combined)
                delete!(terms, word)
            else
                terms[word] = combined
            end
        end
    end
    return PhysicalDysonExpr(terms)
end

function physical_dyson_matrix(a11, a12, a21, a22)
    result = Matrix{PhysicalDysonExpr}(undef, 2, 2)
    result[1, 1] = physical_dyson_expr(a11)
    result[1, 2] = physical_dyson_expr(a12)
    result[2, 1] = physical_dyson_expr(a21)
    result[2, 2] = physical_dyson_expr(a22)
    return result
end

function physical_dyson_matmul(a::Matrix{PhysicalDysonExpr}, b::Matrix{PhysicalDysonExpr})
    result = Matrix{PhysicalDysonExpr}(undef, 2, 2)
    for row in 1:2, column in 1:2
        result[row, column] = a[row, 1] * b[1, column] + a[row, 2] * b[2, column]
    end
    return result
end

function physical_dyson_matsub(a::Matrix{PhysicalDysonExpr}, b::Matrix{PhysicalDysonExpr})
    result = Matrix{PhysicalDysonExpr}(undef, 2, 2)
    for row in 1:2, column in 1:2
        result[row, column] = a[row, column] - b[row, column]
    end
    return result
end

function physical_dyson_scale(factor::Number, matrix::Matrix{PhysicalDysonExpr})
    result = Matrix{PhysicalDysonExpr}(undef, 2, 2)
    for row in 1:2, column in 1:2
        result[row, column] = factor * matrix[row, column]
    end
    return result
end

physical_dyson_identity() = physical_dyson_matrix(1, 0, 0, 1)

function physical_dyson_hs_interaction()
    ψc = physical_dyson_ψ[Classical]
    ψq = physical_dyson_ψ[Quantum]
    χc = physical_dyson_χ[Classical]
    χq = physical_dyson_χ[Quantum]
    forward = ψc^2 * bar(χq) + 2 * ψc * ψq * bar(χc) + (1 // 4) * ψq^2 * bar(χq)
    return ChargedInteractionLagrangian(
        -im * (forward + bar(forward)),
        physical_dyson_ψ => 1,
        physical_dyson_χ => 2;
        parameter=:h,
    )
end

@testset "physical 2PI Dyson closure provenance" begin
    L = physical_dyson_hs_interaction()
    Γ2 = @inferred TwoPIEffectiveAction(L, Val(2), Val(3))
    Σ = @inferred SelfEnergy(Γ2, physical_dyson_ψ)
    Ω = @inferred SelfEnergy(Γ2, physical_dyson_χ)

    @test KC.parameters(Σ) == KC.parameters(Γ2)
    @test KC.parameters(Ω) == KC.parameters(Γ2)
    @test KC.target_family(Σ) == physical_dyson_ψ
    @test KC.target_family(Ω) == physical_dyson_χ

    Σmatrix = @inferred KC.matrix(Σ)
    Ωmatrix = @inferred KC.matrix(Ω)
    @test iszero(Σmatrix[1, 1])
    @test iszero(Ωmatrix[1, 1])
    @test isequal(Σmatrix[1, 2], KC.advanced_component(Σ))
    @test isequal(Σmatrix[2, 1], KC.retarded_component(Σ))
    @test isequal(Σmatrix[2, 2], KC.keldysh_component(Σ))
    @test isequal(Ωmatrix[1, 2], KC.advanced_component(Ω))
    @test isequal(Ωmatrix[2, 1], KC.retarded_component(Ω))
    @test isequal(Ωmatrix[2, 2], KC.keldysh_component(Ω))
end

@testset "atomic left and right Dyson equations" begin
    zero_expr = physical_dyson_scalar(0)
    identity = physical_dyson_identity()

    ℒ = physical_dyson_symbol(:L)
    G_K = physical_dyson_symbol(:G_K)
    G_R = physical_dyson_symbol(:G_R)
    G_A = physical_dyson_symbol(:G_A)
    Σ_A = physical_dyson_symbol(:Sigma_A)
    Σ_R = physical_dyson_symbol(:Sigma_R)
    Σ_K = physical_dyson_symbol(:Sigma_K)

    G = physical_dyson_matrix(G_K, G_R, G_A, zero_expr)
    Σ = physical_dyson_matrix(zero_expr, Σ_A, Σ_R, Σ_K)

    # The source convention is iG₀⁻¹ = [0 L; L 0], hence G₀⁻¹ = -i[0 L; L 0].
    iG0inv = physical_dyson_matrix(zero_expr, ℒ, ℒ, zero_expr)
    G0inv = physical_dyson_scale(-im, iG0inv)
    kernel = physical_dyson_matsub(G0inv, Σ)

    left = physical_dyson_matsub(physical_dyson_matmul(kernel, G), identity)
    @test iszero(left[1, 2])
    @test left[2, 1] == -im * (ℒ * G_K - im * Σ_R * G_K - im * Σ_K * G_A)
    @test left[2, 2] == -im * (ℒ * G_R - im * physical_dyson_scalar(1) - im * Σ_R * G_R)
    @test left[1, 1] == -im * (ℒ * G_A - im * physical_dyson_scalar(1) - im * Σ_A * G_A)

    # The dual Dyson identity keeps the same kernels but reverses all convolution order.
    right = physical_dyson_matsub(physical_dyson_matmul(G, kernel), identity)
    @test iszero(right[2, 1])
    @test right[1, 2] == -im * (G_K * ℒ - im * G_K * Σ_A - im * G_R * Σ_K)
    @test right[1, 1] == -im * (G_R * ℒ - im * physical_dyson_scalar(1) - im * G_R * Σ_R)
    @test right[2, 2] == -im * (G_A * ℒ - im * physical_dyson_scalar(1) - im * G_A * Σ_A)
end

@testset "HS left and right Dyson equations" begin
    g = 3 // 2
    γ = 2 // 5
    λ = complex(g, -γ)
    λbar = conj(λ)
    norm2 = g^2 + γ^2

    zero_expr = physical_dyson_scalar(0)
    one_expr = physical_dyson_scalar(1)
    identity = physical_dyson_identity()

    D_K = physical_dyson_symbol(:D_K)
    D_R = physical_dyson_symbol(:D_R)
    D_A = physical_dyson_symbol(:D_A)
    Ω_A = physical_dyson_symbol(:Omega_A)
    Ω_R = physical_dyson_symbol(:Omega_R)
    Ω_K = physical_dyson_symbol(:Omega_K)

    D = physical_dyson_matrix(D_K, D_R, D_A, zero_expr)
    Ω = physical_dyson_matrix(zero_expr, Ω_A, Ω_R, Ω_K)

    # Corrected Phase-A action kernel. The Dyson inverse is D₀⁻¹ = -i(iD₀⁻¹).
    iD0inv = physical_dyson_matrix(
        zero_expr, -inv(λbar), -inv(λ), complex(0 // 1, -γ / norm2)
    )
    D0inv = physical_dyson_scale(-im, iD0inv)
    kernel = physical_dyson_matsub(D0inv, Ω)

    left = physical_dyson_matsub(physical_dyson_matmul(kernel, D), identity)
    @test iszero(left[1, 2])

    source_D_R = D_R + im * λ * one_expr + im * λ * Ω_R * D_R
    source_D_A = D_A + im * λbar * one_expr + im * λbar * Ω_A * D_A
    source_D_K = D_K + (im * γ / λbar) * D_A + im * λ * Ω_R * D_K + im * λ * Ω_K * D_A

    @test left[2, 2] == (im / λ) * source_D_R
    @test left[1, 1] == (im / λbar) * source_D_A
    @test left[2, 1] == (im / λ) * source_D_K

    # Dual/right Dyson multiplication is an independent operator-order oracle.
    right = physical_dyson_matsub(physical_dyson_matmul(D, kernel), identity)
    @test iszero(right[2, 1])

    source_D_R_right = D_R + im * λ * one_expr + im * λ * D_R * Ω_R
    source_D_A_right = D_A + im * λbar * one_expr + im * λbar * D_A * Ω_A
    source_D_K_right =
        D_K + (im * γ / λ) * D_R + im * λbar * D_K * Ω_A + im * λbar * D_R * Ω_K

    @test right[1, 1] == (im / λ) * source_D_R_right
    @test right[2, 2] == (im / λbar) * source_D_A_right
    @test right[1, 2] == (im / λbar) * source_D_K_right
end
