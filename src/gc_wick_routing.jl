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

@inline function _direct_onepi_wick_pairings(
    args_nc::Vector{Field{S}},
    ::Val{E},
    ::Val{E2};
    regularise=true,
    _set_reg_to_zero=false,
    simplify=false,
) where {S<:Union{Boson,Fermion},E,E2}
    if _use_gc_wick(S, Val(E))
        return _gc_wick_contraction(
            args_nc,
            Val(E),
            Val(E2);
            regularise,
            _set_reg_to_zero,
            simplify,
            onepi_pruning=true,
        )
    end

    pairings = _wick_contraction(
        args_nc, Val(E), Val(E2); regularise, _set_reg_to_zero, simplify
    )
    filter!(entry -> is_irreducible(first(entry).contractions), pairings)
    return pairings
end

function _onepi_wick_contraction!(
    diagrams::Diagrams{C,S,E1,E2},
    a::QMul{A,S};
    regularise=true,
    simplify=false,
    _set_reg_to_zero=false,
) where {C<:Number,A<:Number,S<:Union{Boson,Fermion},E1,E2}
    @assert is_conserved(a)
    @assert is_physical(a)

    pairings = _direct_onepi_wick_pairings(
        a.args_nc, Val(E1), Val(E2); regularise, _set_reg_to_zero, simplify
    )
    make_diagram!(diagrams, pairings, a.arg_c)
    return nothing
end

function _onepi_wick_contraction(
    in_out::QMul{CI,S},
    L::InteractionLagrangian{CL,S},
    ::Val{O},
    ::Val{E};
    simplify=true,
    _set_reg_to_zero=false,
) where {CI<:Number,CL<:Number,S<:Union{Boson,Fermion},O,E}
    @assert number_of_propagators(L) * O + 1 == E "The supplied Val{edges} must equal the interaction's propagator count times Val{order}, plus the external propagator"

    l = length(L.lagrangian)
    C = diagram_coefficient_type(CI, CL)
    diagrams = Diagrams{C,S,E,max_edges(O)}()
    prefactor = wick_prefactor(C, Val(O))
    regularise = should_regularise(L.lagrangian)

    for coefficients in Combinatorics.multiexponents(l, O)
        idxs = indices_from_counts(coefficients)
        mult = Combinatorics.multinomial(coefficients...)
        qmul = mult * prod(L(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))
        term = prefactor * in_out * qmul
        _onepi_wick_contraction!(
            diagrams, term; regularise, simplify, _set_reg_to_zero
        )
    end
    return diagrams
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
