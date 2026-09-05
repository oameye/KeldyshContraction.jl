#################################
#       Wick pairing
#################################

"""A fixed-size Wick pairing with its statistics-dependent exchange sign."""
struct WickPairing{S<:Statistics,E}
    contractions::FixedVector{E,Contraction{S}}
    sign::Int8
end

function WickPairing(
    contractions::Vector{Contraction{S}}, sign::Int8, ::Val{E}
) where {S<:Statistics,E}
    @assert length(contractions) == E "The supplied Val{edges} must match the pairing size"
    fixed = FixedVector{E,Contraction{S}}(c for c in contractions)
    return WickPairing{S,E}(fixed, sign)
end

Base.length(::WickPairing{S,E}) where {S,E} = E
Base.iterate(p::WickPairing) = iterate(p.contractions)
Base.iterate(p::WickPairing, state) = iterate(p.contractions, state)
Base.eltype(::Type{WickPairing{S,E}}) where {S,E} = Contraction{S}

function Diagram(pairing::WickPairing{Boson,E}, ::Val{E}, ::Val{E2}) where {E,E2}
    return Diagram(collect(pairing.contractions), Val(E), Val(E2))
end

pairing_sign(::Type{Boson}, perm) = Int8(1)

"""
    wick_contraction(in_out, L, ::Val{order}, ::Val{edges}; kwargs...)

Compute Wick-contracted diagrams of an interaction with the supplied external fields.
"""
function wick_contraction(
    in_out::QMul, L::InteractionLagrangian, ::Val{O}, ::Val{E}; kwargs...
) where {O,E}
    @assert number_of_propagators(L) * O + 1 == E "The supplied Val{edges} must equal the interaction's propagator count times Val{order}, plus the external propagator"
    return _wick_contraction(in_out, L, Val(E), Val(O); kwargs...)
end

@inline function _wick_contraction(
    in_out::QMul, L::InteractionLagrangian, ::Val{E}, ::Val{O}; kwargs...
) where {E,O}
    l = length(L.lagrangian)

    diagrams = Diagrams{E,max_edges(O)}()
    prefactor = -1 * im * im^O // factorial(O)

    regularise = should_regularise(L.lagrangian)
    for coefficients in Combinatorics.multiexponents(l, O)
        idxs = indices_from_counts(coefficients)
        mult = Combinatorics.multinomial(coefficients...)
        qmul = mult * prod(L(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))
        term = prefactor * in_out * qmul
        wick_contraction!(diagrams, term; regularise, kwargs...)
    end
    return diagrams
end

function wick_contraction(
    in_out::QMul, Ls::LagrangianSum, ::Val{O}, ::Val{E}; simplify::Vector{Bool}, kwargs...
) where {O,E}
    @assert all(number_of_propagators(L) * O + 1 == E for L in arguments(Ls)) "All LagrangianSum terms must produce the supplied number of propagator edges"
    ps = parameters(Ls)
    exponents = Combinatorics.multiexponents(length(Ls), O)
    L_args = arguments(Ls)

    pairs = map(exponents) do coefficients
        idxs = indices_from_counts(coefficients)

        diagrams = if allequal(idxs)
            idx = first(idxs)
            wick_contraction(
                in_out, L_args[idx], Val(O), Val(E); simplify=simplify[idx], kwargs...
            )
        else
            mult = Combinatorics.multinomial(coefficients...)
            qadd = mult * prod(L_args[j](i).lagrangian for (i, j) in enumerate(idxs))
            _simplify = prod(simplify[i] for i in idxs)

            _wick_contraction(in_out, qadd, Val(E), Val(O); simplify=_simplify, kwargs...)
        end
        return (prod(ps[idx] for idx in idxs), diagrams)
    end
    return pairs
end

@inline function _wick_contraction(
    in_out::QMul, a::QAdd, ::Val{E}, ::Val{O}; kwargs...
) where {E,O}
    diagrams = Diagrams{E,max_edges(O)}()
    prefactor = -1 * im * im^O / factorial(O)

    regularise = should_regularise(a)
    for arg in terms(a)
        wick_contraction!(diagrams, prefactor * in_out * arg; regularise, kwargs...)
    end
    return diagrams
end

function wick_contraction!(
    diagrams::Diagrams{E1,E2},
    a::QMul;
    regularise=true,
    simplify=false,
    _set_reg_to_zero=false,
) where {E1,E2}
    @assert is_conserved(a)
    @assert is_physical(a)

    pairings = _wick_contraction(a.args_nc, Val(E1); regularise, _set_reg_to_zero)
    make_diagram!(diagrams, pairings, a.arg_c, simplify)
    return nothing
end

function make_diagram!(
    diagrams::Diagrams{E1,E2},
    pairings::Vector{WickPairing{Boson,E1}},
    arg_c,
    simplify::Bool,
) where {E1,E2}
    isempty(pairings) && return nothing
    imag_factor = im^E1
    for pairing in pairings
        diagram, prefactor = make_diagram_pair(
            pairing, arg_c, imag_factor, simplify, Val(E1), Val(E2)
        )
        push!(diagrams, diagram, prefactor)
    end
    return nothing
end

function make_diagram_pair(
    pairing::WickPairing{Boson,E}, arg_c, imag_factor, simplify::Bool, ::Val{E}, ::Val{E2}
) where {E,E2}
    contractions = collect(pairing.contractions)
    contractions′, prefactor =
        simplify ? advanced_to_retarded(contractions, arg_c) : (contractions, arg_c)
    prefactor *= pairing.sign
    sort!(contractions′; by=sort_by_position_and_type)
    edges = FixedVector{E,Edge}(Edge(contraction) for contraction in contractions′)
    return Diagram(edges, Val(E2)) => imag_factor * prefactor
end

"""
    prepare_args(args, ::Val{E})

Split a canonical field product into `E` unbarred fields and the reversed sequence of
`E` barred fields. Wick permutations act on indices of that reversed partner sequence.
This ordering is the parity reference used by statistics-dependent `pairing_sign` methods.
"""
function prepare_args(args::Vector{Field{S}}, ::Val{E}) where {S<:Statistics,E}
    @assert length(args) == 2E "Number of fields must be twice the pairing size"
    destroys = args[1:E]
    creates = reverse(args[(E + 1):end])
    return destroys, creates
end

function _wick_contraction(
    args_nc::Vector{Field{S}}, ::Val{E}; regularise=true, _set_reg_to_zero=false
)::Vector{WickPairing{S,E}} where {S<:Statistics,E}
    destroys, creates = prepare_args(args_nc, Val(E))
    ps = map(position, args_nc)
    skip = has_in(ps) && has_out(ps)

    wick_pairings = WickPairing{S,E}[]

    for perm in SmallCombinatorics.permutations(E)
        if skip && isone(first(perm))
            continue
        end
        contractions, fail = wick_contract(
            destroys, creates, perm; regularise, _set_reg_to_zero
        )

        if fail || !is_connected(contractions) || has_zero_loop(contractions)
            continue
        end

        canonical = canonicalize(contractions)
        push!(wick_pairings, WickPairing(canonical, pairing_sign(S, perm), Val(E)))
    end

    return wick_pairings
end

function wick_contract(
    destroys::Vector{Field{S}},
    creates::Vector{Field{S}},
    perm;
    regularise=true,
    _set_reg_to_zero=false,
) where {S<:Statistics}
    contractions = Contraction{S}[]
    fail = false
    for (k, l) in pairs(perm)
        potential = Contraction(destroys[k], creates[l])
        if !contraction_filter(potential)
            fail = true
            break
        end
        if regularise
            if !regular(potential)
                fail = true
                break
            end
            different_position = !allequal(position.(potential))
            if _set_reg_to_zero && (different_position || is_keldysh(potential))
                potential = map(set_reg_to_zero, potential)
            end
        end
        push!(contractions, potential)
    end
    return contractions, fail
end

######################
# Vacuum Contractions
######################

function _wick_contraction(a::QAdd, ::Val{E}; kwargs...) where {E}
    args = terms(a)
    @assert all(number_of_propagators(arg) == E for arg in args)
    @assert is_bulk(a) "The private two-argument _wick_contraction entry point is for vacuum terms"

    diagrams = Diagrams{E,topology_length(E + 1)}()
    regularise = should_regularise(a)
    for arg in args
        wick_contraction!(diagrams, arg; regularise, kwargs...)
    end
    return diagrams
end

function _wick_contraction(a::QMul, ::Val{E}; kwargs...) where {E}
    @assert is_conserved(a)
    @assert is_physical(a)
    @assert number_of_propagators(a) == E
    @assert is_bulk(a) "The private two-argument _wick_contraction entry point is for vacuum terms"

    diagrams = Diagrams{E,topology_length(E + 1)}()
    regularise = should_regularise(a)
    wick_contraction!(diagrams, a; regularise, kwargs...)
    return diagrams
end
