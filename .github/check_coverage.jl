function source_relative_path(path::AbstractString)
    normalized = normpath(path)
    root = normpath(pwd())
    if startswith(normalized, root)
        return relpath(normalized, root)
    end
    marker = string(Base.Filesystem.path_separator, "src", Base.Filesystem.path_separator)
    index = findfirst(marker, normalized)
    isnothing(index) && return normalized
    return normalized[(first(index) + 1):end]
end

function uncovered_source_lines(path::AbstractString)
    uncovered = Dict{String,Vector{Int}}()
    source = nothing

    for line in eachline(path)
        if startswith(line, "SF:")
            source = source_relative_path(line[4:end])
        elseif startswith(line, "DA:") && source !== nothing && startswith(source, "src/")
            fields = split(line[4:end], ',')
            line_number = parse(Int, fields[1])
            count = parse(Int, fields[2])
            iszero(count) || continue
            push!(get!(uncovered, source, Int[]), line_number)
        end
    end

    return uncovered
end

coverage_file = length(ARGS) == 1 ? only(ARGS) : "lcov.info"
uncovered = uncovered_source_lines(coverage_file)

if isempty(uncovered)
    println("Source line coverage: 100%")
    exit(0)
end

println(stderr, "Source line coverage is below 100%. Uncovered package-owned lines:")
for source in sort!(collect(keys(uncovered)))
    lines = sort!(unique(uncovered[source]))
    println(stderr, "  ", source, ": ", join(lines, ", "))
end
exit(1)
