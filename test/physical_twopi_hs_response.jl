using KeldyshContraction, Test
using SymbolicUtils
import KeldyshContraction as KC

@qfields hs_response_ψ::Boson hs_response_χ::Boson

function hs_response_interaction()
    ψc = hs_response_ψ[Classical]
    ψq = hs_response_ψ[Quantum]
    χc = hs_response_χ[Classical]
    χq = hs_response_χ[Quantum]
    forward = ψc^2 * bar(χq) + 2 * ψc * ψq * bar(χc) + (1 // 4) * ψq^2 * bar(χq)
    return ChargedInteractionLagrangian(
        -im * (forward + bar(forward)), hs_response_ψ => 1, hs_response_χ => 2; parameter=:h
    )
end

hs_response_R(λ, ΩR) = (; num=(-(λ^2) * ΩR), den=1 + λ * ΩR)
hs_response_A(λbar, ΩA) = (; num=(-(λbar^2) * ΩA), den=1 + λbar * ΩA)

function hs_response_F(λ, λbar, γ, ΩR, ΩA, ΩF)
    den = (1 + λ * ΩR) * (1 + λbar * ΩA)
    num = -(λ * λbar) * ΩF + γ * (λ * ΩR + λbar * ΩA + λ * λbar * ΩR * ΩA)
    return (; num, den)
end

function hs_response_F_compact(λ, λbar, γ, ΩR, ΩA, ΩF)
    den = (1 + λ * ΩR) * (1 + λbar * ΩA)
    num = γ * den - γ - λ * λbar * ΩF
    return (; num, den)
end

function hs_response_iszero(expression)
    expanded = SymbolicUtils.expand(expression)
    simplified = SymbolicUtils.simplify(expanded)
    return SymbolicUtils._iszero(simplified)
end

@testset "homogeneous regular HS Dyson response" begin
    @syms λ λbar γ ΩR ΩA ΩF

    DR = hs_response_R(λ, ΩR)
    DA = hs_response_A(λbar, ΩA)
    DF = hs_response_F(λ, λbar, γ, ΩR, ΩA, ΩF)

    # Retarded/advanced Kadanoff--Baym equations, with denominators cleared exactly.
    @test hs_response_iszero(DR.num + λ^2 * ΩR * DR.den + λ * ΩR * DR.num)
    @test hs_response_iszero(DA.num + λbar^2 * ΩA * DA.den + λbar * ΩA * DA.num)

    # The statistical equation is checked before introducing any commuting rational-function
    # production type. At homogeneous Fourier order all factors are scalar; clearing the
    # common denominator leaves an exact polynomial identity.
    QR = DR.den
    QA = DA.den
    @test isequal(DF.den, QR * QA)
    statistical_residual =
        λbar * DF.num - (
            -γ * DA.num * QR - λ * λbar^2 * ΩF * DF.den + γ * λ * λbar * ΩR * DF.den -
            λ * λbar * ΩR * DF.num - λ * λbar * ΩF * DA.num * QR
        )
    @test hs_response_iszero(statistical_residual)

    # The compact spelling is algebraically identical, not a second response ansatz.
    compact = hs_response_F_compact(λ, λbar, γ, ΩR, ΩA, ΩF)
    @test isequal(compact.den, DF.den)
    @test hs_response_iszero(compact.num - DF.num)

    # The regular response contains no remnant of the singular D0 contact.
    λsample0 = complex(3 // 2, -2 // 5)
    λbarsample0 = conj(λsample0)
    @test hs_response_R(λsample0, 0) == (; num=0, den=1)
    @test hs_response_A(λbarsample0, 0) == (; num=0, den=1)
    @test hs_response_F(λsample0, λbarsample0, 2 // 5, 0, 0, 0) == (; num=0, den=1)

    # Elastic limit γ -> 0.
    @syms g
    elastic = hs_response_F(g, g, 0, ΩR, ΩA, ΩF)
    @test hs_response_iszero(elastic.num + g^2 * ΩF)
    @test hs_response_iszero(elastic.den - (1 + g * ΩR) * (1 + g * ΩA))

    # Physical R/A conjugation at an exact complex sample point.
    λsample = complex(7 // 5, -3 // 8)
    ΩRsample = complex(2 // 7, 5 // 11)
    ΩAsample = conj(ΩRsample)
    retarded = hs_response_R(λsample, ΩRsample)
    advanced = hs_response_A(conj(λsample), ΩAsample)
    @test advanced.num == conj(retarded.num)
    @test advanced.den == conj(retarded.den)
end

@testset "HS response is tied to the χ-target 2PI polarization" begin
    Γ2 = @inferred TwoPIEffectiveAction(hs_response_interaction(), Val(2), Val(3))
    Ω = @inferred SelfEnergy(Γ2, hs_response_χ)
    ΩF = @inferred KC.twopi_fourier_self_energy(Γ2, hs_response_χ)

    @test target_family(Ω) == hs_response_χ
    @test target_family(ΩF) == hs_response_χ
    @test parameters(Ω) == parameters(Γ2)
    @test parameters(ΩF) == parameters(Γ2)

    # Cutting the χ line from the physical G²D skeleton leaves the pair bubble: two ψ lines,
    # one external pair momentum, and one independent loop momentum. C2 therefore consumes
    # the compiler-generated Ω rather than introducing a separate polarization formula.
    for component in
        (KC.keldysh_component(ΩF), KC.retarded_component(ΩF), KC.advanced_component(ΩF))
        @test !isempty(component)
        for (graph, contributions) in component
            @test graph.external_count == 1
            @test graph.loop_count == 1
            families = [
                KC.field_family(edge.out) for edge in KC.contractions(graph.coordinate)
            ]
            @test length(families) == 2
            @test all(==(hs_response_ψ), families)
            @test !isempty(contributions)
        end
    end
end
