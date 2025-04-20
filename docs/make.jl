CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using KeldyshContraction
using Documenter

using Plots
default(; fmt=:png)
# Gotta set this environment variable when using the GR run-time on CI machines.
# This happens as examples will use Plots.jl to make plots and movies.
# See: https://github.com/jheinen/GR.jl/issues/278
ENV["GKSwstype"] = "100"

include("pages.jl")

makedocs(;
    sitename="KeldyshContraction.jl",
    authors="Orjan Ameyep",
    modules=KeldyshContraction,
    format=Documenter.HTML(;
        canonical="https://oameye.github.io/KeldyshContraction.jl"
    ),
    pages=pages,
    clean=true,
    linkcheck=true,
    warnonly=:missing_docs,
    draft=false,#,(!CI),
    doctest=false,  # We test it in the CI, no need to run it here
    checkdocs=:exports,
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
