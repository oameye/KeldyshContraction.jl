using Graphs

g = DiGraph(Edge.([(2, 1), (2,3), (2,3), (3, 2), (1, 3)]))

g = DiGraph(Edge.([(1, 2), (1,3), (2,3)]))
edges(g) |> collect

is_strongly_connected(g)

simplecycles(g)
simplecycles_hawick_james(g)
circuit(g)
