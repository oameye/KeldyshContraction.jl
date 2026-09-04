#################################
#       Contraction
#################################

"""
    wick_contraction(in_out::QMul, L::InteractionLagrangian, ::Val{order}, ::Val{edges}; kwargs...)
    wick_contraction(in_out::QMul, Ls::LagrangianSum, ::Val{order}, ::Val{edges}; simplify, kwargs...)

Compute all possible Wick contractions of the interaction Lagrangian with the external
fields in `in_out`.

Wick contractions decompose products of quantum field operators into sums of products
of propagators (two-point correlation functions). The rules of the contraction are:
  - Conservation (equal numbers of creation/annihilation operators)
  - Physicality (proper time ordering)
  - No quantum-quantum contractions
  - If the fields have a [`Regularisation`](@ref) applied, the contractions are
    regularised.

The function returns a new expression of propagators as `Diagrams`. For a
`LagrangianSum`, one parameter/diagram pair is returned for each parameter monomial.
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
        # TODO: remove complex conjugate to go from 10 to only 6 terms
        idxs = indices_from_counts(coefficients) # will be of length order
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
        idxs = indices_from_counts(coefficients) # will be of length order

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
    for arg in arguments(a)
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

    contractions = _wick_contraction(a.args_nc; regularise, _set_reg_to_zero)
    make_diagram!(diagrams, contractions, a.arg_c, simplify)
    return nothing
end

function make_diagram!(
    diagrams::Diagrams{E1,E2}, contractions, arg_c, simplify::Bool
) where {E1,E2}
    if isempty(contractions)
        return nothing
    end
    number_of_contractions = length(first(contractions))
    imag_factor = im^(number_of_contractions) # Contraction becomes propagator
    # foreach(contractions) do c
    for c in contractions
        diagram, prefactor = make_diagram_pair(
            c, arg_c, imag_factor, simplify, Val(E1), Val(E2)
        )
        push!(diagrams, diagram, prefactor)
    end
    return nothing
end

function make_diagram_pair(
    c, arg_c, imag_factor, simplify::Bool, ::Val{E}, ::Val{E2}
) where {E,E2}
    c′, prefactor = simplify ? advanced_to_retarded(c, arg_c) : (c, arg_c)
    sort!(c′; by=sort_by_position_and_type)
    edges = FixedVector{E,Edge}(Edge(contraction) for contraction in c′)
    return Diagram(edges, Val(E2)) => imag_factor * prefactor
end

"""
We split up the fields into two groups, `destroys` and `creates`. We can can combute all
possible pairs by permutating the create vector. To avoid pairing up the In() and Out()
fields, we have made sure that the destroy and create vectors are ordered with the in and
out fields first. Computing the permutatins in lexicographic order, we can skip the first
(n-1)! permutations.
"""
function _wick_contraction(
    args_nc::Vector{<:QField}; regularise=true, _set_reg_to_zero=false
)::Vector{Vector{Contraction}}
    destroys, creates, n_destroy = prepare_args(args_nc)

    number_of_combinations = factorial(n_destroy)
    ps = map(position, args_nc)
    skip = has_in(ps) && has_out(ps)

    wick_contractions = Vector{Contraction}[]

    for (i, perm) in enumerate(SmallCombinatorics.permutations(n_destroy))
        if skip && isone(first(perm))
            continue # skip the in-out contraction
        end
        contraction, fail = _wick_contract(
            destroys, creates, perm; regularise, _set_reg_to_zero
        )

        # TODO ∨ You can probably cache this
        if fail || !is_connected(contraction) || has_zero_loop(contraction)
            continue
        else
            canonical_ordered_contraction = canonicalize(contraction)
            push!(wick_contractions, canonical_ordered_contraction)
        end
    end

    return wick_contractions
end
function _wick_contract(destroys, creates, perm; regularise=true, _set_reg_to_zero=false)
    contraction = Contraction[]
    fail = false
    for (k, l) in pairs(perm)
        potential_contraction = (destroys[k], creates[l])
        if !contraction_filter(potential_contraction)
            fail = true
            break
        end
        if regularise
            if !regular(potential_contraction)
                fail = true
                break
            else
                same_position = !allequal(position.(potential_contraction))
                if _set_reg_to_zero && (same_position || is_keldysh(potential_contraction))
                    potential_contraction = set_reg_to_zero.(potential_contraction)
                end
            end
        end
        push!(contraction, potential_contraction)
    end
    return contraction, fail
end
function prepare_args(args::Vector{<:QField})
    _length = length(args)
    @assert _length % 2 == 0 "Number of fields must be even"

    n_destroy = _length ÷ 2

    destroys = args[1:n_destroy]
    creates = reverse(args[(n_destroy + 1):end])
    return destroys, creates, n_destroy
end

######################
# Vacuum Contractions
######################

function _wick_contraction(a::QAdd, ::Val{E}; kwargs...) where {E}
    args = SymbolicUtils.arguments(a)
    @assert all(number_of_propagators(arg) == E for arg in args)
    if is_bulk(a) # for vacuum calculations
        diagrams = Diagrams{E,topology_length(E + 1)}()
    else
        @warn """
        The private `_wick_contraction` helper is intended for vacuum calculations only.
        Instead, use `wick_contraction` with InteractionLagrangian.
        """
        diagrams = Diagrams{E,topology_length(E)}()
    end

    regularise = should_regularise(a)
    for arg in args
        wick_contraction!(diagrams, arg; regularise, kwargs...)
    end
    return diagrams
end # keep for vacuum calculations
function _wick_contraction(a::QMul, ::Val{E}; kwargs...) where {E}
    @assert is_conserved(a)
    @assert is_physical(a)

    @assert number_of_propagators(a) == E
    if is_bulk(a) # for vacuum calculations
        diagrams = Diagrams{E,topology_length(E + 1)}()
    else
        @warn """
        The private `_wick_contraction` helper is intended for vacuum calculations only.
        Instead, use `wick_contraction` with InteractionLagrangian.
        """
        diagrams = Diagrams{E,topology_length(E)}()
    end

    regularise = should_regularise(a)
    wick_contraction!(diagrams, a; regularise, kwargs...)
    return diagrams
end # keep for vacuum calculations

# The following were used to check for bugs, we leave them here for reference
# but they are not used in the main code.
# function check_sorted(args)
#     args′ = sort(args; by=position)
#     args′′ = sort(args′; by=ladder)
#     @assert isequal(args, args′′) "Arguments are not sorted"
# end
# function check_to_many_bulk(contraction, args_nc)
#     pos = map(x -> Int(position(x)), Iterators.flatten(contraction))
#     for i in unique(pos)
#         if count(x -> x == i, pos) > 4
#             @show args_nc
#             error("Contraction is not unique")
#         end
#     end
# end
