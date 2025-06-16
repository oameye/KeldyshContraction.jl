
#################################
#       Contraction
#################################
"""
    wick_contraction(expr::QTerm)

Compute all possible Wick contractions of quantum fields in the expression `expr`.

Wick contractions decompose products of quantum field operators into sums of products
of propagators (two-point correlation functions). The rules of the contraction are:
  - Conservation (equal numbers of creation/annihilation operators)
  - Physicality (proper time ordering)
  - No quantum-quantum contractions
  - If the fields have a [`Regularisation`](@ref) applied, the contractions are
    regularised. The [`Regularisation`](@ref) property is set to zero after the reguralisation.

The function handles two types of inputs:
- `QAdd`: Distributes the contraction over sums
- `QMul`: Contracts products of fields into propagators

The function returns a new expression of propagators of type `SymbolicUtils.Symbol`.

"""
function wick_contraction(in_out::QMul, L::InteractionLagrangian, order::Int64; kwargs...)
    l = length(L.lagrangian)

    E = number_of_propagators(L) * order + 1 # +1 for in_out
    diagrams = Diagrams{E,topology_length(E)}()
    prefactor = -1 * im * im^order / factorial(order)

    regularise = should_regularise(L.lagrangian)
    for coefficients in Combinatorics.multiexponents(l, order)
        # TODO: remove complex conjugate to go from 10 to only 6 terms
        idxs = indices_from_counts(coefficients) # will be of length order
        mult = Combinatorics.multinomial(coefficients...)
        qmul = mult * prod(L(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))
        wick_contraction!(diagrams, prefactor * in_out * qmul; regularise, kwargs...)
    end
    return diagrams
end

function wick_contraction!(diagrams::Diagrams, a::QMul; regularise=true, simplify=true)
    @assert is_conserved(a)
    @assert is_physical(a)

    contractions = _wick_contraction(a.args_nc; regularise)
    make_diagram!(diagrams, contractions, a.arg_c, simplify)
    return nothing
end
function make_diagram!(diagrams::Diagrams, contractions, arg_c, simplify::Bool)
    if isempty(contractions)
        return nothing
    end
    number_of_contractions = length(first(contractions))
    imag_factor = im^(number_of_contractions) # Contraction becomes propagator
    # foreach(contractions) do c
    for c in contractions
        diagram, prefactor = make_diagram_pair(c, arg_c, imag_factor, simplify)
        push!(diagrams, diagram, prefactor)
    end
    return nothing
end
function make_diagram_pair(c, arg_c, imag_factor, simplify::Bool)
    c′, prefactor = simplify ? advanced_to_retarded(c, arg_c) : (c, arg_c)
    return Diagram(c′) => imag_factor * prefactor
end

"""
We split up the fields into two groups, `destroys` and `creates`. We can can combute all
possible pairs by permutating the create vector. To avoid pairing up the In() and Out()
fields, we have made sure that the destroy and create vectors are ordered with the in and
out fields first. Computing the permutatins in lexicographic order, we can skip the first
(n-1)! permutations.
"""
function _wick_contraction(
    args_nc::Vector{<:QField}; regularise=true
)::Vector{Vector{Contraction}}
    destroys, creates, n_destroy = prepare_args(args_nc)

    number_of_combinations = factorial(n_destroy)
    ps = map(position, args_nc)
    skip = has_in(ps) && has_out(ps)

    wick_contractions = Vector{Contraction}[]

    for (i, perm) in enumerate(SmallCollections.Combinatorics.permutations(n_destroy))
        if skip && isone(first(perm))
            continue # skip the in-out contraction
        end
        contraction, fail = _wick_contract(destroys, creates, perm; regularise)

        # TODO ∨ You can probably cache this
        if fail || !is_connected(contraction) || has_zero_loop(contraction)
            continue
        else
            push!(wick_contractions, canonicalize(contraction))
        end
    end
    return wick_contractions
end
function _wick_contract(destroys, creates, perm; regularise=true)
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
                potential_contraction = set_reg_to_zero.(potential_contraction)
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

function _wick_contraction(a::QAdd; kwargs...)::Diagrams
    args = SymbolicUtils.arguments(a)
    E = number_of_propagators(first(args))
    if is_bulk(a) # for vacuum calculations
        diagrams = Diagrams{E,topology_length(E + 1)}()
    else # known type unstable
        @warn """
        Directly using `wick_contraction` on a `QAdd` is only publicly supported for vacuum
        calculations. Instead, use `wick_contraction` with InteractionLagrangian.
        """
        diagrams = Diagrams{E,topology_length(E)}()
    end

    regularise = should_regularise(a)
    for arg in args
        wick_contraction!(diagrams, arg; regularise, kwargs...)
    end
    return diagrams
end # keep for vacuum calculations
function _wick_contraction(a::QMul; kwargs...)::Diagrams
    @assert is_conserved(a)
    @assert is_physical(a)

    E = number_of_propagators(a)
    if is_bulk(a) # for vacuum calculations
        diagrams = Diagrams{E,topology_length(E + 1)}()
    else # known type unstable
        @warn """
        Directly using `wick_contraction` on a `QAdd` is only publicly supported for vacuum
        calculations. Instead, use `wick_contraction` with InteractionLagrangian.
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
