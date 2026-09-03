using Combinatorics, BenchmarkTools
import KeldyshContraction as KC

perm = [1, 2, 3, 4, 5]
KC.next_permutation!(perm, 1, 5)
perm
function f(n)
    number_of_combinations = factorial(n)
    to_skip = factorial(n - 1) # due in-out contraction constraint
    iter = 1:n
    contraction = Vector{Int}[]
    for i in (to_skip + 1):number_of_combinations
        perm = Combinatorics.nthperm(iter, i)
        push!(contraction, perm)
    end
    return contraction
end

function g(n)
    number_of_combinations = factorial(n)
    to_skip = factorial(n - 1) # due in-out contraction constraint
    iter = 1:n
    return contraction = collect(permutations(iter))[(to_skip + 1):number_of_combinations]
end

function h(n)
    number_of_combinations = factorial(n)
    to_skip = factorial(n - 1) # due in-out contraction constraint
    iter = 1:n
    contraction = Vector{Int}[]

    perm = Combinatorics.nthperm(iter, to_skip)
    for i in (to_skip + 1):number_of_combinations
        KC.next_permutation!(perm, 1, n)
        push!(contraction, collect(perm))
    end
    return contraction
end

f(5)
g(5)
h(5)

n = 9
@btime f($n)
@btime g($n)
@btime h($n)
# ^ f is faster than g
