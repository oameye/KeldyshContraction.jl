CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using KeldyshContraction
using Documenter

using Plots: Plots
Plots.default(; fmt=:png)
# Gotta set this environment variable when using the GR run-time on CI machines.
# This happens as examples will use Plots.jl to make plots and movies.
# See: https://github.com/jheinen/GR.jl/issues/278
ENV["GKSwstype"] = "100"

include("pages.jl")
include("make_md_examples.jl")

# The README.md file is used as the deployed documentation home page. Avoid copying it during
# local LiveServer builds, where doing so would create a rebuild loop.
if CI
    cp(
        normpath(@__FILE__, "../../README.md"),
        normpath(@__FILE__, "../src/index.md");
        force=true,
    )
end

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
