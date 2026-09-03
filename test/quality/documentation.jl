using KeldyshContraction
using Test
using Documenter

@testset "Documentation" begin
    DocMeta.setdocmeta!(
        KeldyshContraction, :DocTestSetup, :(using KeldyshContraction); recursive=true
    )
    Documenter.doctest(KeldyshContraction)
end
