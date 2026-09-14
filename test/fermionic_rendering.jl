using KeldyshContraction, Test

@qfields render_ψ::Fermion
ψ₁, ψ₂ = render_ψ[One], render_ψ[Two]

@testset "fermionic LaTeX rendering" begin
    @test repr(MIME"text/latex"(), ψ₁) == "\$render_ψ^1\$"
    @test repr(MIME"text/latex"(), ψ₂) == "\$render_ψ^2\$"
    @test repr(MIME"text/latex"(), bar(ψ₁)) == "\$\\bar{render_ψ^1}\$"
    @test repr(MIME"text/latex"(), bar(ψ₂)) == "\$\\bar{render_ψ^2}\$"

    expression = ψ₁ * ψ₂ * bar(ψ₂) * bar(ψ₁)
    latex = repr(MIME"text/latex"(), expression)
    @test contains(latex, "render_ψ^1")
    @test contains(latex, "render_ψ^2")
    @test contains(latex, "\\bar{render_ψ^1}")
    @test contains(latex, "\\bar{render_ψ^2}")
end
