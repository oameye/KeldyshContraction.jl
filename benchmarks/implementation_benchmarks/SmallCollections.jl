using BenchmarkTools

using SmallCombinatorics: SmallCombinatorics
using Combinatorics: Combinatorics
using Combinat: Combinat
n = 9;
@btime sum(p[1] for p in SmallCombinatorics.permutations($n)); # SmallCollections.jl
# 13.424 ms (10 allocations: 192 bytes)
@btime sum(p[1] for p in Combinatorics.permutations(1:($n))); # Combinatorics.jl
# 5.469 s (725761 allocations: 55.37 MiB)
@btime sum(p[1] for p in Combinat.Permutations($n)); # Combinat.jl
# 7.113 ms (362881 allocations: 44.30 MiB)
