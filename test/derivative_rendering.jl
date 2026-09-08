using KeldyshContraction, Test
using KeldyshContraction: write_latex

@qfields render_ϕ::Boson render_ψ::Fermion

@testset "derivative field rendering" begin
    c = render_ϕ[Classical]
    ψ₁ = render_ψ[One]

    ∂xyc = partial(partial(c, :y), :x)
    ∂xψ = partial(ψ₁, :x)

    @test repr(∂xyc) == "∂x∂yrender_ϕᶜ"
    @test repr(bar(∂xyc)) == "∂x∂ȳrender_ϕᶜ"
    @test repr(∂xψ) == "∂xrender_ψ"

    @test repr(MIME"text/latex"(), ∂xyc) == "\$\\partial_{x} \\partial_{y} render_ϕ^c\$"
    @test repr(MIME"text/latex"(), bar(∂xyc)) ==
        "\$\\partial_{x} \\partial_{y} \\bar{render_ϕᶜ}\$"

    io = IOBuffer()
    @test @inferred(show(io, ∂xψ)) === nothing
    io = IOBuffer()
    @test @inferred(write_latex(io, ∂xyc)) === nothing
end
