using KeldyshContraction, BenchmarkTools
using Combinatorics
using KeldyshContraction: In, Out, Classical, Quantum, Plus, Minus, arguments
import KeldyshContraction as KC

@qfields ϕ::Boson
ϕᶜ, ϕᴾ = ϕ[Classical], ϕ[Quantum]

L_int =
    im * (
        0.5 * bar(ϕᶜ) * bar(ϕᴾ) * (ϕᶜ(Minus) * ϕᶜ(Minus) + ϕᴾ(Minus) * ϕᴾ(Minus)) -
        0.5 * ϕᶜ(Plus) * ϕᴾ(Plus) * (bar(ϕᶜ) * bar(ϕᶜ) + bar(ϕᴾ) * bar(ϕᴾ)) +
        bar(ϕᶜ) * bar(ϕᴾ) * (ϕᶜ(Plus) * ϕᴾ(Plus) + ϕᶜ(Minus) * ϕᴾ(Minus))
    )

expr = ϕᶜ(Out()) * bar(ϕᶜ)(In()) * L_int

function _wick_contraction(
    args_nc::Vector{Field{Boson}}
)::Vector{Vector{Vector{Field{Boson}}}}
    _partitions = Combinatorics.partitions(args_nc, length(args_nc) ÷ 2)
    wick_contractions = Vector{Vector{Field{Boson}}}[]
    for v in _partitions
        if _contraction_filter(v)
            push!(wick_contractions, v)
        end
    end
    return wick_contractions
end

function _contraction_filter(v)
    istwo = all(length.(v) .== 2) # only two-point contractions
    if !istwo
        return false
    elseif !all(KC.is_conserved.(v)) # contractions with barred/unbarred fields
        return false
    elseif !all(KC.is_physical_propagator.(v)) # propagators are physical
        return false
    else # there should be no quantum-quantum contractions
        return !KC.has_qq_contraction(v)
    end
end

to_bench = arguments(expr)[1].args_nc
@btime wick_contraction($to_bench) # 38.851 μs (369 allocations: 15.67 KiB)
@btime _wick_contraction($to_bench) # 57.920 μs (1608 allocations: 82.83 KiB)
