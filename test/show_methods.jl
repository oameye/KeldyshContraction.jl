using KeldyshContraction, Test
using KeldyshContraction: In, Out, Edge, Bulk
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus

@qfields ϕ::Boson(Classical) ψ::Boson(Quantum)
@qfields ϕᶜ::Boson(Classical) ϕᴾ::Boson(Quantum)

@testset "Symbols" begin
    input = [
        ϕ,
        bar(ϕ),
        bar(ϕ(Plus)),
        bar(ϕ(Minus)),
        Edge(ϕ, bar(ϕ)(In())),
        Edge(ϕ(Out()), bar(ϕ)),
        Edge(ϕ(Out()), bar(ψ)),
        Edge(ψ(Out()), bar(ϕ)),
        Edge(ψ(Bulk(2)), bar(ϕ)),
    ]
    output = [
        "ϕ",
        "̄ϕ",
        "̄ϕ⁺",
        "̄ϕ⁻",
        "Gᴷ(y₁,x₂)",
        "Gᴷ(x₁,y₁)",
        "Gᴿ(x₁,y₁)",
        "Gᴬ(x₁,y₁)",
        "Gᴬ(y₂,y₁)",
    ]
    for (i, o) in zip(input, output)
        @test sprint(show, i) == o
        @test repr(i) == o
    end
    s = IOBuffer(sizehint=0)
    @inferred show(s, ϕ * bar(ϕ))

    output_latex = [
        "\$\\phi\$",
        "\$\\bar{ϕ}\$",
        "\$\\bar{ϕ}^+\$",
        "\$\\bar{ϕ}^{-}\$",
        "\$G^K\\left( y_1, x_2 \\right)\$",
        "\$G^K\\left( x_1, y_1 \\right)\$",
        "\$G^R\\left( x_1, y_1 \\right)\$",
        "\$G^A\\left( x_1, y_1 \\right)\$",
    ]

    for (i, o) in zip(input, output_latex)
        @test sprint(show, MIME"text/latex"(), i) == o
        @test repr(MIME"text/latex"(), i) == o
    end
end

@testset "Term" begin
    @qfields ϕᶜ::Boson(Classical) ϕᴾ::Boson(Quantum)

    L_int = im * (0.5 * bar(ϕᶜ) * bar(ϕᴾ) * (ϕᶜ * ϕᶜ))
    @test repr(L_int) == "(0.0 + 0.5im)*(ϕᶜ*ϕᶜ*̄ϕᴾ*̄ϕᶜ)"
    @test repr(MIME"text/latex"(), L_int) ==
        "\$0.5 i \\phi^c \\phi^c \\bar{\\phi^P} \\bar{\\phi^c}\$"
end

@testset "Structs" begin
    using KeldyshContraction: Diagram, Diagrams
    using SymbolicUtils
    @syms g::Number

    L = InteractionLagrangian(ϕ * ψ * bar(ϕ) * bar(ψ))
    @test repr(MIME"text/plain"(), L) ==
        "Interaction Lagrangian with fields ϕ and ψ:\n(ψ*ϕ*̄ψ*̄ϕ)"

    @test repr(MIME"text/latex"(), L) == "\$\\psi \\phi \\bar{\\psi} \\bar{\\phi}\$"

    ds = Diagrams(
        [Diagram([Edge(ϕ, bar(ϕ))], Val(1), Val(0))], Complex{Rational{Int}}(1.0)
    )
    DP = DressedPropagator(ds, ds, ds, 1, g)
    @test repr(MIME"text/plain"(), DP) ==
        "Dressed Propagator:\nkeldysh:  Gᴷ(y₁,y₁)\nretarded: Gᴷ(y₁,y₁)\nadvanced: Gᴷ(y₁,y₁)"
end

@testset "Momentum" begin
    using KeldyshContraction: Momenta, Edge, Bulk, Momentum

    @test repr(Momenta(0)) == "k"
    @test repr(Momenta(1)) == "q₁"
    @test repr(Momenta(2)) == "q₂"

    ms = Momenta([1, 1, -1], [Momentum(1), Momentum(2), Momentum(0)])
    @test repr(ms) == "q₁ + q₂ - k"

    e = Edge(ψ(Bulk(2)), bar(ϕ))
    e0 = Edge(e, Momenta(0))
    e1 = Edge(e, Momenta(1))
    e2 = Edge(e, Momenta(2))

    @test repr(e0) == "Gᴬ(k)"
    @test repr(e1) == "Gᴬ(q₁)"
    @test repr(e2) == "Gᴬ(q₂)"
end
