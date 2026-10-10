# Production crossover for the symmetry-aware GraphCombinations Wick backend.
#
# The direct matcher remains the low-order path. Production-shaped certification shows mixed
# results at five edges: GC wins for the coherent bosonic workload, but the colored two-body-loss
# workload allocates more, and low-order fermions likewise do not improve robustly. From seven
# edges onward GC wins in both time and allocations across coherent, lossy, plain-fermion, and
# derivative-fermion workloads.

@inline _use_gc_wick(::Type{Boson}, ::Val{E}) where {E} = E >= 7
@inline _use_gc_wick(::Type{Fermion}, ::Val{E}) where {E} = E >= 7

@inline function _production_wick_pairings(
    args_nc::Vector{Field{S}},
    ::Val{E},
    ::Val{E2};
    regularise=true,
    _set_reg_to_zero=false,
    simplify=false,
) where {S<:Union{Boson,Fermion},E,E2}
    if _use_gc_wick(S, Val(E))
        return _gc_wick_contraction(
            args_nc, Val(E), Val(E2); regularise, _set_reg_to_zero, simplify
        )
    end
    return _wick_contraction(
        args_nc, Val(E), Val(E2); regularise, _set_reg_to_zero, simplify
    )
end

function wick_contraction!(
    diagrams::Diagrams{C,S,E1,E2},
    a::QMul{A,S};
    regularise=true,
    simplify=false,
    _set_reg_to_zero=false,
) where {C<:Number,A<:Number,S<:Union{Boson,Fermion},E1,E2}
    @assert is_conserved(a)
    @assert is_physical(a)

    pairings = _production_wick_pairings(
        a.args_nc, Val(E1), Val(E2); regularise, _set_reg_to_zero, simplify
    )
    make_diagram!(diagrams, pairings, a.arg_c)
    return nothing
end
