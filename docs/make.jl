CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using KeldyshContraction
using Documenter

include("pages.jl")
include("make_md_examples.jl")

makedocs(;
    sitename="KeldyshContraction.jl",
    authors="Orjan Ameye",
    modules=KeldyshContraction,
    format=Documenter.HTML(; canonical="https://oameye.github.io/KeldyshContraction.jl"),
    pages=pages,
    pagesonly=true,
    clean=true,
    linkcheck=true,
    draft=false,
    doctest=false,  # Doctests are run in the test suite.
    checkdocs=:public,
)

if CI
    deploydocs(;
        repo="github.com/oameye/KeldyshContraction.jl",
        devbranch="main",
        target="build",
        branch="gh-pages",
        push_preview=true,
    )
end
