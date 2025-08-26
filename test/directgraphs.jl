using Graphs

g = DiGraph(Edge.([(1, 2), (3,2), (3,2), (2, 3), (3, 1)]))

edges(g) |> collect

simplecycles(g)

c = Channel{Vector{Int}}()
task = @async Graphs.itercycles(g, c)
while !istaskdone(task)
    cycle = take!(c)
    # if length(cycle) == nv(g)
    #     println(cycle)
    # end
    yield()
end

simplecycles_hawick_james(g)
