using KeldyshContraction, Test
using KeldyshContraction: In, Out, Edge, Bulk
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus

@qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)

@qfields ϕᶜ::Destroy(Classical) ϕᴾ::Destroy(Quantum)

@testset "Symbols" begin
    input = [
        ϕ,
        ϕ',
        ϕ(Plus)',
        ϕ(Minus)',
        Edge(ϕ, ϕ'(In())),
        Edge(ϕ(Out()), ϕ'),
        Edge(ϕ(Out()), ψ'),
        Edge(ψ(Out()), ϕ'),
        Edge(ψ(Bulk(2)), ϕ'),
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
    @inferred show(s, ϕ*ϕ')

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
    @qfields ϕᶜ::Destroy(Classical) ϕᴾ::Destroy(Quantum)

    L_int = im*(0.5 * ϕᶜ' * ϕᴾ' * (ϕᶜ * ϕᶜ))
    @test repr(L_int) == "(0.0 + 0.5im)*(ϕᶜ*ϕᶜ*̄ϕᴾ*̄ϕᶜ)"
    @test repr(MIME"text/latex"(), L_int) ==
        "\$0.5 i \\phi^c \\phi^c \\bar{\\phi^P} \\bar{\\phi^c}\$"
end

@testset "Structs" begin
    using KeldyshContraction: Diagram, Diagrams
    using SymbolicUtils

    L = InteractionLagrangian(ϕ * ψ * ϕ' * ψ')
    @test repr(MIME"text/plain"(), L) ==
        "Interaction Lagrangian with fields ϕ and ψ:\n(ϕ*ψ*̄ψ*̄ϕ)"

    @test repr(MIME"text/latex"(), L) == "\$\\phi \\psi \\bar{\\psi} \\bar{\\phi}\$"

    ds = Diagrams([Diagram([Edge(ϕ, ϕ')])], complex(1.0))
    DP = DressedPropagator(ds, ds, ds, 1, SymbolicUtils.Sym{Number}(:g))
    @test repr(MIME"text/plain"(), DP) ==
        "Dressed Propagator:\nkeldysh:  Gᴷ(y₁,y₁)\nretarded: Gᴷ(y₁,y₁)\nadvanced: Gᴷ(y₁,y₁)"
end

@testset "Momentum" begin
    using KeldyshContraction: Momenta, Edge, Bulk
    e = Edge(ψ(Bulk(2)), ϕ')
    e0 = Edge(e, Momenta(0))
    e1 = Edge(e, Momenta(1))
    e2 = Edge(e, Momenta(2))

    @test repr(e0) == "Gᴬ(k)"
    @test repr(e1) == "Gᴬ(q₁)"
    @test repr(e2) == "Gᴬ(q₂)"
end
