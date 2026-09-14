using Literate

const EXAMPLES_IN = joinpath(@__DIR__, "..", "examples")
const OUTPUT_MD_DIR = joinpath(@__DIR__, "src", "examples")
const EXAMPLES = [
    joinpath(EXAMPLES_IN, "bosonic_transport.jl"),
    joinpath(EXAMPLES_IN, "fermionic_transport.jl"),
]

mkpath(OUTPUT_MD_DIR)

if isempty(get(ENV, "CI", ""))
    # Needed only for local builds; CI supplies repository URLs automatically.
    extra_literate_config = Dict(
        "repo_root_path" => abspath(joinpath(@__DIR__, "..")),
        "repo_root_url" => "file://" * abspath(joinpath(@__DIR__, "..")),
    )
else
    extra_literate_config = Dict()
end

for example in EXAMPLES
    Literate.markdown(
        example,
        OUTPUT_MD_DIR;
        flavor=Literate.DocumenterFlavor(),
        config=extra_literate_config,
    )
end
