using KeldyshContraction, Test
using KeldyshContraction: In, Out, Edge, Bulk, write_latex
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus

@qfields ϕ::Boson
ϕᶜ, ϕᴾ = ϕ[Classical], ϕ[Quantum]

@testset "Symbols" begin
    input = [
        ϕᶜ,
        bar(ϕᶜ),
        bar(ϕᶜ(Plus)),
        bar(ϕᶜ(Minus)),
        Edge(ϕᶜ, bar(ϕᶜ)(In())),
        Edge(ϕᶜ(Out()), bar(ϕᶜ)),
        Edge(ϕᶜ(Out()), bar(ϕᴾ)),
        Edge(ϕᴾ(Out()), bar(ϕᶜ)),
        Edge(ϕᴾ(Bulk(2)), bar(ϕᶜ)),
    ]
    output = [
        "ϕᶜ",
        "̄ϕᶜ",
        "̄ϕᶜ⁺",
        "̄ϕᶜ⁻",
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
    @test @inferred(show(s, ϕᶜ * bar(ϕᶜ))) === nothing

    output_latex = [
        "\$\\phi^c\$",
        "\$\\bar{ϕᶜ}\$",
        "\$\\bar{ϕᶜ}^+\$",
        "\$\\bar{ϕᶜ}^{-}\$",
        "\$G^K\\left( y_1, x_2 \\right)\$",
        "\$G^K\\left( x_1, y_1 \\right)\$",
        "\$G^R\\left( x_1, y_1 \\right)\$",
        "\$G^A\\left( x_1, y_1 \\right)\$",
        "\$G^A\\left( y_2, y_1 \\right)\$",
    ]

    for (i, o) in zip(input, output_latex)
        @test sprint(show, MIME"text/latex"(), i) == o
        @test repr(MIME"text/latex"(), i) == o
    end

    @test contains(repr(MIME"text/latex"(), Edge(ϕᶜ(Plus), bar(ϕᶜ))), "^+")
    @test contains(repr(MIME"text/latex"(), Edge(ϕᶜ(Minus), bar(ϕᶜ))), "^-")

    io = IOBuffer()
    @test @inferred(show(io, MIME"text/latex"(), ϕᶜ)) === nothing
    io = IOBuffer()
    @test @inferred(write_latex(io, ϕᶜ)) === nothing
end

@testset "Term" begin
    L_int = im * (0.5 * bar(ϕᶜ) * bar(ϕᴾ) * (ϕᶜ * ϕᶜ))
    @test repr(L_int) == "(0.0 + 0.5im)*(ϕᶜ*ϕᶜ*̄ϕᴾ*̄ϕᶜ)"
    @test repr(MIME"text/latex"(), L_int) ==
        "\$0.5 i \\phi^c \\phi^c \\bar{\\phi^P} \\bar{\\phi^c}\$"

    negative_term = -(ϕᶜ * bar(ϕᶜ))
    negative_term_latex = repr(MIME"text/latex"(), negative_term)
    @test startswith(negative_term_latex, "\$- ")

    difference = ϕᶜ * bar(ϕᶜ) - ϕᴾ * bar(ϕᴾ)
    difference_latex = repr(MIME"text/latex"(), difference)
    @test contains(difference_latex, " - ")
    @test !contains(difference_latex, "+ -")

    io = IOBuffer()
    @test @inferred(show(io, L_int)) === nothing
    io = IOBuffer()
    @test @inferred(show(io, MIME"text/latex"(), L_int)) === nothing
    io = IOBuffer()
    @test @inferred(write_latex(io, L_int)) === nothing
    io = IOBuffer()
    @test @inferred(write_latex(io, difference)) === nothing
end

@testset "Public LaTeX coefficient rendering" begin
    using SymbolicUtils
    @qfields a::Boson
    @syms g::Number
    aᶜ = a[Classical]
    @test repr(MIME"text/latex"(), aᶜ) == "\$a^c\$"

    fields = ϕᶜ * bar(ϕᶜ)
    @test contains(repr(MIME"text/latex"(), (1 // 2) * fields), "\\frac{1}{2}")
    @test contains(repr(MIME"text/latex"(), im * fields), "i ")
    @test contains(repr(MIME"text/latex"(), -im * fields), "-i ")
    @test contains(repr(MIME"text/latex"(), 2im * fields), "2 i ")
    @test contains(repr(MIME"text/latex"(), complex(1.0, 2.0) * fields), "(1.0 + 2.0 i)")
    @test contains(repr(MIME"text/latex"(), complex(1.0, -2.0) * fields), "(1.0 - 2.0 i)")

    symbolic_sum = g * fields + ϕᴾ * bar(ϕᴾ)
    symbolic_latex = repr(MIME"text/latex"(), symbolic_sum)
    @test contains(symbolic_latex, "g")
    @test contains(symbolic_latex, " + ")
end

@testset "Structs" begin
    using KeldyshContraction: Diagram, Diagrams
    using SymbolicUtils
    @syms g::Number

    L = InteractionLagrangian(ϕᶜ * ϕᴾ * bar(ϕᶜ) * bar(ϕᴾ))
    @test repr(MIME"text/plain"(), L) ==
        "Interaction Lagrangian with fields ϕ:\n(ϕᴾ*ϕᶜ*̄ϕᴾ*̄ϕᶜ)"

    @test repr(MIME"text/latex"(), L) == "\$\\phi^P \\phi^c \\bar{\\phi^P} \\bar{\\phi^c}\$"

    io = IOBuffer()
    @test @inferred(show(io, L)) === nothing
    io = IOBuffer()
    @test @inferred(show(io, MIME"text/latex"(), L)) === nothing
    io = IOBuffer()
    @test @inferred(write_latex(io, L)) === nothing

    ds = Diagrams(
        [Diagram([Edge(ϕᶜ, bar(ϕᶜ))], Val(1), Val(0))], Complex{Rational{Int}}(1.0)
    )
    @test repr(MIME"text/latex"(), ds) == "\$G^K\\left( y_1, y_1 \\right)\$"

    diagram = first(keys(ds.diagrams))
    scaled_ds = Diagrams([diagram], Complex{Rational{Int}}(2))
    scaled_latex = repr(MIME"text/latex"(), scaled_ds)
    @test contains(scaled_latex, "2 G^K")
    @test contains(repr(scaled_ds), "Gᴷ")

    negative_ds = Diagrams([diagram], Complex{Rational{Int}}(-2))
    negative_latex = repr(MIME"text/latex"(), negative_ds)
    @test contains(negative_latex, "- 2 G^K")
    @test !contains(negative_latex, "+ -")

    two_edge_diagram = Diagram([Edge(ϕᶜ, bar(ϕᶜ)), Edge(ϕᶜ, bar(ϕᶜ))], Val(2), Val(0))
    two_edge_latex = repr(MIME"text/latex"(), two_edge_diagram)
    @test count(==('G'), two_edge_latex) == 2

    DP = DressedPropagator(ds, ds, ds, Val(1), parameter_monomial(g))
    @test repr(MIME"text/plain"(), DP) ==
        "Dressed Propagator:\nkeldysh:  Gᴷ(y₁,y₁)\nretarded: Gᴷ(y₁,y₁)\nadvanced: Gᴷ(y₁,y₁)"

    io = IOBuffer()
    @test @inferred(show(io, diagram)) === nothing
    io = IOBuffer()
    @test @inferred(show(io, MIME"text/latex"(), diagram)) === nothing
    io = IOBuffer()
    @test @inferred(write_latex(io, diagram)) === nothing
end

@testset "Momentum" begin
    using KeldyshContraction: Momenta, Edge, Bulk, Momentum

    @test repr(Momenta(0)) == "k"
    @test repr(Momenta(1)) == "q₁"
    @test repr(Momenta(2)) == "q₂"

    ms = Momenta([1, 1, -1], [Momentum(1), Momentum(2), Momentum(0)])
    @test repr(ms) == "q₁ + q₂ - k"

    zero_ms = Momenta([0], [Momentum(0)])
    negative_ms = Momenta([-1], [Momentum(1)])

    e = Edge(ϕᴾ(Bulk(2)), bar(ϕᶜ))
    e0 = Edge(e, Momenta(0))
    e1 = Edge(e, Momenta(1))
    e2 = Edge(e, Momenta(2))
    em = Edge(e, ms)
    ezero = Edge(e, zero_ms)
    enegative = Edge(e, negative_ms)

    @test repr(e0) == "Gᴬ(k)"
    @test repr(e1) == "Gᴬ(q₁)"
    @test repr(e2) == "Gᴬ(q₂)"

    @test repr(MIME"text/latex"(), e0) == "\$G^A\\left( k \\right)\$"
    @test repr(MIME"text/latex"(), e1) == "\$G^A\\left( q_1 \\right)\$"
    @test repr(MIME"text/latex"(), e2) == "\$G^A\\left( q_2 \\right)\$"
    @test repr(MIME"text/latex"(), em) == "\$G^A\\left( q_1 + q_2 - k \\right)\$"
    @test repr(MIME"text/latex"(), ezero) == "\$G^A\\left( 0 \\right)\$"
    @test repr(MIME"text/latex"(), enegative) == "\$G^A\\left( -q_1 \\right)\$"

    io = IOBuffer()
    @test @inferred(show(io, ms)) === nothing
    io = IOBuffer()
    @test @inferred(write_latex(io, ms)) === nothing
end
