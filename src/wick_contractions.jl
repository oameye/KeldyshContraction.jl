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

function Diagram(pairing::WickPairing{S,E}, ::Val{E}, ::Val{E2}) where {S<:Statistics,E,E2}
    contractions = Contraction{S}[contraction for contraction in pairing.contractions]
    return Diagram(contractions, Val(E), Val(E2))
end

pairing_sign(::Type{Boson}, perm) = Int8(1)

function wick_prefactor(::Type{C}, ::Val{O}) where {C<:Number,O}
    exact = convert(ComplexRationals, -im * im^O) / factorial(O)
    return convert(C, exact)
end

"""
    wick_contraction(in_out, L, ::Val{order}, ::Val{edges}; kwargs...)

Compute Wick-contracted diagrams of an interaction with the supplied external fields.
"""
function wick_contraction(
    in_out::QMul{CI,S}, L::InteractionLagrangian{CL,S}, ::Val{O}, ::Val{E}; kwargs...
) where {CI<:Number,CL<:Number,S<:Statistics,O,E}
    @assert number_of_propagators(L) * O + 1 == E "The supplied Val{edges} must equal the interaction's propagator count times Val{order}, plus the external propagator"
    return _wick_contraction(in_out, L, Val(E), Val(O); kwargs...)
end

@inline function _wick_contraction(
    in_out::QMul{CI,S}, L::InteractionLagrangian{CL,S}, ::Val{E}, ::Val{O}; kwargs...
) where {CI<:Number,CL<:Number,S<:Statistics,E,O}
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
        wick_contraction!(diagrams, term; regularise, kwargs...)
    end
    return diagrams
end

function wick_contraction(
    in_out::QMul{CI,S},
    Ls::LagrangianSum{CL,S},
    ::Val{O},
    ::Val{E};
    simplify::Vector{Bool},
    kwargs...,
) where {CI<:Number,CL<:Number,S<:Statistics,O,E}
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
    in_out::QMul{CI,S}, a::QAdd{CA,S}, ::Val{E}, ::Val{O}; kwargs...
) where {CI<:Number,CA<:Number,S<:Statistics,E,O}
    C = diagram_coefficient_type(CI, CA)
    diagrams = Diagrams{C,S,E,max_edges(O)}()
    prefactor = wick_prefactor(C, Val(O))

    regularise = should_regularise(a)
    for arg in terms(a)
        wick_contraction!(diagrams, prefactor * in_out * arg; regularise, kwargs...)
    end
    return diagrams
end

function wick_contraction!(
    diagrams::Diagrams{C,S,E1,E2},
    a::QMul{A,S};
    regularise=true,
    simplify=false,
    _set_reg_to_zero=false,
) where {C<:Number,A<:Number,S<:Statistics,E1,E2}
    @assert is_conserved(a)
    @assert is_physical(a)

    pairings = _wick_contraction(
        a.args_nc, Val(E1), Val(E2); regularise, _set_reg_to_zero, simplify
    )
    make_diagram!(diagrams, pairings, a.arg_c)
    return nothing
end

function make_diagram!(
    diagrams::Diagrams{C,S,E1,E2},
    pairings::Vector{Tuple{WickPairing{S,E1},FixedVector{E2,Int},Int}},
    arg_c,
) where {C<:Number,S<:Statistics,E1,E2}
    isempty(pairings) && return nothing
    imag_factor = convert(C, im^E1)
    for (pairing, topology, multiplicity) in pairings
        diagram, prefactor = make_diagram_pair(
            pairing, topology, multiplicity, arg_c, imag_factor, Val(E1), Val(E2)
        )
        push!(diagrams, diagram, prefactor)
    end
    return nothing
end

function make_diagram_pair(
    pairing::WickPairing{S,E},
    topology::FixedVector{E2,Int},
    multiplicity::Int,
    arg_c,
    imag_factor,
    ::Val{E},
    ::Val{E2},
) where {S<:Statistics,E,E2}
    contractions = Contraction{S}[contraction for contraction in pairing.contractions]
    sort!(contractions; by=sort_by_position_and_type)
    edges = FixedVector{E,Edge{S}}(Edge(contraction) for contraction in contractions)
    prefactor = arg_c * pairing.sign * multiplicity
    return Diagram{S,E,E2}(edges, topology) => imag_factor * prefactor
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

const WICK_ORBIT_EXHAUSTIVE_LIMIT = 5

@inline function wick_contraction_isless(a::Contraction{S}, b::Contraction{S}) where {S}
    isequal(a.out, b.out) || return isless(a.out, b.out)
    return isless(a.in, b.in)
end

function wick_fixed_isless(
    a::FixedVector{E,Contraction{S}}, b::FixedVector{E,Contraction{S}}
) where {S<:Statistics,E}
    for i in 1:E
        isequal(a[i], b[i]) && continue
        return wick_contraction_isless(a[i], b[i])
    end
    return false
end

function sorted_wick_key(
    contractions::Vector{Contraction{S}}, ::Val{E}
) where {S<:Statistics,E}
    sort!(contractions; lt=wick_contraction_isless)
    return FixedVector{E,Contraction{S}}(contractions)
end

const WickFieldKey = Tuple{Int,Int8,Int,Symbol,Int,NTuple{4,Int16}}
const WickContractionKey = Tuple{WickFieldKey,WickFieldKey}

@inline function wick_field_key(f::Field, mapping::Dict{Position,Position})::WickFieldKey
    p = position(f)
    mapped_position = is_bulk(p) ? index(mapping[p]) : index(p)
    return (
        Int(orientation(f)),
        mapped_position,
        Int(keldysh_index(f)),
        name(f),
        Int(regularisation(f)),
        slots(field_indices(f)),
    )
end

@inline function wick_compact_contraction_key(
    c::Contraction, mapping::Dict{Position,Position}
)::WickContractionKey
    return (wick_field_key(c.out, mapping), wick_field_key(c.in, mapping))
end

function wick_compact_isless(
    a::Vector{WickContractionKey}, b::Vector{WickContractionKey}
)::Bool
    @inbounds for i in eachindex(a)
        isequal(a[i], b[i]) && continue
        return isless(a[i], b[i])
    end
    return false
end

function foreach_position_relabeling!(
    f::F,
    positions::Vector{Position},
    mapping::Dict{Position,Position},
    first_label::Int,
    k::Int,
) where {F}
    if k > length(positions)
        f(mapping)
        return nothing
    end

    for j in k:length(positions)
        positions[k], positions[j] = positions[j], positions[k]
        old_position = positions[k]
        mapping[old_position] = Bulk(first_label + k - 1)
        foreach_position_relabeling!(f, positions, mapping, first_label, k + 1)
        delete!(mapping, old_position)
        positions[k], positions[j] = positions[j], positions[k]
    end
    return nothing
end

function foreach_wick_relabeling!(
    f::F, anchors::Vector{Position}, others::Vector{Position}
) where {F}
    mapping = Dict{Position,Position}()
    foreach_position_relabeling!(anchors, mapping, 1, 1) do anchor_mapping
        foreach_position_relabeling!(
            others, anchor_mapping, length(anchors) + 1, 1
        ) do complete_mapping
            return f(complete_mapping)
        end
    end
    return nothing
end

wick_orbit_key(::Type{S}, raw::FixedVector{E,Contraction{S}}) where {S<:Statistics,E} = raw

function wick_orbit_key(::Type{Boson}, raw::FixedVector{E,Contraction{Boson}}) where {E}
    contractions = Contraction{Boson}[contraction for contraction in raw]
    graph_positions = canonicalization_positions(contractions)
    bulk_positions = Position[position for position in graph_positions if is_bulk(position)]

    if length(bulk_positions) > WICK_ORBIT_EXHAUSTIVE_LIMIT
        return sorted_wick_key(contractions, Val(E))
    end

    anchors = out_bulk_positions(contractions)
    sort!(anchors)
    others = Position[position for position in bulk_positions if position ∉ anchors]

    scratch = Vector{WickContractionKey}(undef, E)
    best = similar(scratch)
    best_mapping = Dict{Position,Position}()
    best_set = false
    foreach_wick_relabeling!(anchors, others) do mapping
        @inbounds for i in 1:E
            scratch[i] = wick_compact_contraction_key(contractions[i], mapping)
        end
        sort!(scratch)
        if !best_set || wick_compact_isless(scratch, best)
            copyto!(best, scratch)
            empty!(best_mapping)
            for (old_position, new_position) in mapping
                best_mapping[old_position] = new_position
            end
            best_set = true
        end
    end

    relabeled = Contraction{Boson}[
        relabel_bulk_positions(contraction, best_mapping) for contraction in contractions
    ]
    return sorted_wick_key(relabeled, Val(E))
end

"""
Build the locally admissible partners for each destroying field.

Family/index compatibility, orientation, physicality, QQ exclusion, and regularisation are
pair-local constraints. Applying them once here avoids constructing all `E!` permutations
only to reject almost all of them later.
"""
function wick_candidates(
    destroys::Vector{Field{S}},
    creates::Vector{Field{S}},
    ::Val{E};
    regularise=true,
    _set_reg_to_zero=false,
    skip_external_pair=false,
) where {S<:Statistics,E}
    return ntuple(Val(E)) do k
        candidates = Tuple{Int,Contraction{S}}[]
        for l in 1:E
            skip_external_pair && k == 1 && l == 1 && continue

            potential = Contraction(destroys[k], creates[l])
            contraction_filter(potential) || continue
            regularise && !regular(potential) && continue

            if regularise
                different_position = !allequal(position.(potential))
                if _set_reg_to_zero && (different_position || is_keldysh(potential))
                    potential = map(set_reg_to_zero, potential)
                end
            end
            push!(candidates, (l, potential))
        end
        return candidates
    end
end

function _foreach_wick_matching!(
    f::F,
    candidates::Tuple,
    contractions::Vector{C},
    permutation::Vector{Int},
    used::Vector{Bool},
    k::Int,
) where {F,C<:Contraction}
    if k > length(candidates)
        f(contractions, permutation)
        return nothing
    end

    for (l, contraction) in candidates[k]
        used[l] && continue
        used[l] = true
        contractions[k] = contraction
        permutation[k] = l
        _foreach_wick_matching!(f, candidates, contractions, permutation, used, k + 1)
        used[l] = false
    end
    return nothing
end

function foreach_wick_matching(
    f::F, candidates::Candidates, ::Val{E}
) where {F,E,Candidates<:Tuple}
    candidate_vector_type = fieldtype(Candidates, 1)
    candidate_type = fieldtype(eltype(candidate_vector_type), 2)
    contractions = Vector{candidate_type}(undef, E)
    permutation = Vector{Int}(undef, E)
    used = fill(false, E)
    _foreach_wick_matching!(f, candidates, contractions, permutation, used, 1)
    return nothing
end

function wick_matching_key(
    contractions::Vector{C}, candidates::Candidates, ::Val{E}
)::NTuple{E,UInt8} where {C<:Contraction,Candidates<:Tuple,E}
    return ntuple(Val(E)) do k
        idx = findfirst(candidates[k]) do candidate
            return isequal(candidate[2], contractions[k])
        end
        isnothing(idx) &&
            error("Wick matching contraction is not present in its candidate row")
        return UInt8(idx)
    end
end

function contractions_from_matching_key(
    key::NTuple{E,UInt8}, candidates::Candidates, ::Type{Contraction{S}}
) where {S<:Statistics,E,Candidates<:Tuple}
    contractions = Vector{Contraction{S}}(undef, E)
    @inbounds for k in 1:E
        contractions[k] = candidates[k][Int(key[k])][2]
    end
    return contractions
end

function _wick_contraction(
    args_nc::Vector{Field{S}}, ::Val{E}; regularise=true, _set_reg_to_zero=false
)::Vector{WickPairing{S,E}} where {S<:Statistics,E}
    destroys, creates = prepare_args(args_nc, Val(E))
    ps = map(position, args_nc)
    skip = has_in(ps) && has_out(ps)
    candidates = wick_candidates(
        destroys, creates, Val(E); regularise, _set_reg_to_zero, skip_external_pair=skip
    )

    wick_pairings = WickPairing{S,E}[]
    foreach_wick_matching(candidates, Val(E)) do contractions, permutation
        passes_wick_filters(contractions) || return nothing

        canonical = canonicalize(contractions)
        push!(wick_pairings, WickPairing(canonical, pairing_sign(S, permutation), Val(E)))
        return nothing
    end
    return wick_pairings
end

"""
Generate final physical Wick-pairing representatives together with their static legacy topology.

Matching permutations are first merged by a compact exact raw-contraction signature, including the
statistics-dependent signed multiplicity. The signature indexes the precomputed candidate rows and
therefore avoids storing full fixed contraction tuples as dictionary keys. Connectivity, causal
filtering, optional advanced-to-retarded simplification, and color-aware physical canonicalization
are then paid once per unique raw pairing. Legacy uncolored topology is computed once per final
physical diagram.
"""
function _wick_contraction(
    args_nc::Vector{Field{S}},
    ::Val{E},
    ::Val{E2};
    regularise=true,
    _set_reg_to_zero=false,
    simplify=false,
)::Vector{Tuple{WickPairing{S,E},FixedVector{E2,Int},Int}} where {S<:Statistics,E,E2}
    destroys, creates = prepare_args(args_nc, Val(E))
    ps = map(position, args_nc)
    skip = has_in(ps) && has_out(ps)
    candidates = wick_candidates(
        destroys, creates, Val(E); regularise, _set_reg_to_zero, skip_external_pair=skip
    )

    matching_weights = Dict{NTuple{E,UInt8},Int}()
    foreach_wick_matching(candidates, Val(E)) do contractions, permutation
        key = wick_matching_key(contractions, candidates, Val(E))
        weight = Int(pairing_sign(S, permutation))
        matching_weights[key] = get(matching_weights, key, 0) + weight
        return nothing
    end

    canonical_weights = Dict{FixedVector{E,Contraction{S}},Int}()
    for (key, weight) in matching_weights
        iszero(weight) && continue
        contractions = contractions_from_matching_key(key, candidates, Contraction{S})
        passes_wick_filters(contractions) || continue

        final_contractions, simplification_sign = if simplify
            advanced_to_retarded(contractions, 1)
        else
            contractions, 1
        end
        canonical = canonicalize(final_contractions)
        canonical_key = sorted_wick_key(canonical, Val(E))
        final_weight = weight * Int(simplification_sign)
        canonical_weights[canonical_key] =
            get(canonical_weights, canonical_key, 0) + final_weight
    end

    wick_pairings = Tuple{WickPairing{S,E},FixedVector{E2,Int},Int}[]
    sizehint!(wick_pairings, length(canonical_weights))
    for (canonical, weight) in canonical_weights
        iszero(weight) && continue
        contractions = Contraction{S}[contraction for contraction in canonical]
        topology = legacy_topology(contractions, Val(E2))
        pairing = WickPairing(contractions, Int8(sign(weight)), Val(E))
        push!(wick_pairings, (pairing, topology, abs(weight)))
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
# Explicit-shape contractions
######################

function _wick_contraction(
    a::QAdd{C,S}, ::Val{E}, ::Val{E2}; kwargs...
) where {C<:Number,S<:Statistics,E,E2}
    args = terms(a)
    @assert all(number_of_propagators(arg) == E for arg in args)

    D = diagram_coefficient_type(C)
    diagrams = Diagrams{D,S,E,E2}()
    regularise = should_regularise(a)
    for arg in args
        wick_contraction!(diagrams, arg; regularise, kwargs...)
    end
    return diagrams
end

function _wick_contraction(
    a::QMul{C,S}, ::Val{E}, ::Val{E2}; kwargs...
) where {C<:Number,S<:Statistics,E,E2}
    @assert is_conserved(a)
    @assert is_physical(a)
    @assert number_of_propagators(a) == E

    D = diagram_coefficient_type(C)
    diagrams = Diagrams{D,S,E,E2}()
    regularise = should_regularise(a)
    wick_contraction!(diagrams, a; regularise, kwargs...)
    return diagrams
end

######################
# Vacuum Contractions
######################

function _wick_contraction(
    a::QAdd{C,S}, ::Val{E}; kwargs...
) where {C<:Number,S<:Statistics,E}
    @assert is_bulk(a) "The private two-argument _wick_contraction entry point is for vacuum terms"
    return _wick_contraction(a, Val(E), Val(topology_length(E + 1)); kwargs...)
end

function _wick_contraction(
    a::QMul{C,S}, ::Val{E}; kwargs...
) where {C<:Number,S<:Statistics,E}
    @assert is_bulk(a) "The private two-argument _wick_contraction entry point is for vacuum terms"
    return _wick_contraction(a, Val(E), Val(topology_length(E + 1)); kwargs...)
end
