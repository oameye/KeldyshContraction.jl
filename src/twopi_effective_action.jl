"""Formal full-propagator vacuum diagram used inside the 2PI effective action."""
struct TwoPIDiagram{S<:Statistics,E,E2}
    contractions::FixedVector{E,Contraction{S}}
    topology::FixedVector{E2,Int}
end

function Base.isequal(a::TwoPIDiagram{S,E,E2}, b::TwoPIDiagram{S,E,E2}) where {S,E,E2}
    return isequal(a.contractions, b.contractions)
end
Base.:(==)(a::TwoPIDiagram{S,E,E2}, b::TwoPIDiagram{S,E,E2}) where {S,E,E2} = isequal(a, b)
Base.hash(d::TwoPIDiagram, h::UInt) = hash(TwoPIDiagram, hash(d.contractions, h))

twopi_contractions(d::TwoPIDiagram) = d.contractions
twopi_topology(d::TwoPIDiagram) = d.topology

"""
$(DocStringExtensions.TYPEDEF)

Fixed-order two-particle-irreducible vacuum contribution ``Γ₂`` generated from a Keldysh
interaction action.

Unlike ordinary perturbative propagator diagrams, this representation is a functional of formal
full propagators. It therefore retains Keldysh components that vanish only after restricting to the
physical propagator manifold, including the quantum--quantum component. Physical causal pruning is
not applied before functional differentiation.

`TwoPIEffectiveAction` is diagrammatic compiler IR. It does not solve a Dyson or Kadanoff--Baym
equation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoPIEffectiveAction{C<:Number,S<:Statistics,O,E,E2}
    "Canonical connected 2PI vacuum diagrams with exact coefficients."
    diagrams::Dict{TwoPIDiagram{S,E,E2},C}
    "Canonical perturbation-parameter monomial."
    parameter::ParameterMonomial
    "Physical field families appearing in the interaction."
    families::Vector{FieldFamily{S}}
end

order(::TwoPIEffectiveAction{C,S,O}) where {C,S,O} = O
statistics(::TwoPIEffectiveAction{C,S}) where {C,S} = S
parameters(Γ::TwoPIEffectiveAction) = Γ.parameter
field_families(Γ::TwoPIEffectiveAction) = copy(Γ.families)
Base.length(Γ::TwoPIEffectiveAction) = length(Γ.diagrams)
Base.iszero(Γ::TwoPIEffectiveAction) = isempty(Γ.diagrams)
Base.iterate(Γ::TwoPIEffectiveAction) = iterate(Γ.diagrams)
Base.iterate(Γ::TwoPIEffectiveAction, state) = iterate(Γ.diagrams, state)

"""Return a copy of the exact canonical ``Γ₂`` diagram map."""
twopi_terms(Γ::TwoPIEffectiveAction) = copy(Γ.diagrams)

function _twopi_formal_candidates(
    destroys::Vector{Field{S}}, creates::Vector{Field{S}}, ::Val{E}
) where {S<:Statistics,E}
    return ntuple(Val(E)) do k
        candidates = Tuple{Int,Contraction{S}}[]
        for l in 1:E
            potential = Contraction(destroys[k], creates[l])
            contraction_compatible(S, potential.out, potential.in) || continue
            balanced_orientation(potential) || continue
            is_physical_propagator(potential) || continue
            push!(candidates, (l, potential))
        end
        return candidates
    end
end

function _twopi_formal_pairings(
    args_nc::Vector{Field{S}}, ::Val{E}, ::Val{E2}
) where {S<:Statistics,E,E2}
    destroys, creates = prepare_args(args_nc, Val(E))
    candidates = _twopi_formal_candidates(destroys, creates, Val(E))

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
        is_two_particle_irreducible(contractions) || continue
        canonical = canonicalize(contractions)
        canonical_key = sorted_wick_key(canonical, Val(E))
        canonical_weights[canonical_key] = get(canonical_weights, canonical_key, 0) + weight
    end

    result = Tuple{TwoPIDiagram{S,E,E2},Int}[]
    sizehint!(result, length(canonical_weights))
    for (canonical, weight) in canonical_weights
        iszero(weight) && continue
        contractions = Contraction{S}[contraction for contraction in canonical]
        topology = legacy_topology(contractions, Val(E2))
        diagram = TwoPIDiagram{S,E,E2}(canonical, topology)
        push!(result, (diagram, weight))
    end
    return result
end

function _twopi_uniform_field_count(L::InteractionLagrangian)
    counts = unique(length(term.args_nc) for term in terms(L.lagrangian))
    length(counts) == 1 || throw(
        ArgumentError(
            "TwoPIEffectiveAction currently requires interaction monomials with a common field count",
        ),
    )
    return only(counts)
end

"""
    TwoPIEffectiveAction(L::InteractionLagrangian, ::Val{order}, ::Val{edges})

Generate the connected 2PI vacuum contribution ``Γ₂`` at fixed interaction order directly from
`L`. `edges` is the number of full propagator lines in each vacuum diagram, so for a uniform
`N`-field interaction it must satisfy `2edges == order*N`.

The generator uses formal full-propagator Wick matching. In particular, it does not impose
`G_qq = 0`, causal vacuum-loop cancellation, equal-time retarded/advanced simplification, or any
other physical-propagator reduction before the effective action is differentiated.
"""
function TwoPIEffectiveAction(
    L::InteractionLagrangian{C,S}, ::Val{O}, ::Val{E}
) where {C<:Number,S<:Statistics,O,E}
    O > 0 || throw(ArgumentError("2PI interaction order must be positive"))
    field_count = _twopi_uniform_field_count(L)
    field_count * O == 2E || throw(
        ArgumentError(
            "2PI vacuum edge count must satisfy 2*edges == order*fields_per_vertex"
        ),
    )

    D = diagram_coefficient_type(C)
    E2 = max_edges(O)
    diagrams = Dict{TwoPIDiagram{S,E,E2},D}()
    prefactor = wick_prefactor(D, Val(O))
    nterms = length(L.lagrangian)

    for coefficients in Combinatorics.multiexponents(nterms, O)
        idxs = indices_from_counts(coefficients)
        multiplicity = Combinatorics.multinomial(coefficients...)
        vertex_product =
            multiplicity * prod(L(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))
        term = prefactor * vertex_product

        for (diagram, wick_weight) in _twopi_formal_pairings(term.args_nc, Val(E), Val(E2))
            value = _simplify(convert(D, im^E) * convert(D, term.arg_c) * wick_weight)
            combined = _simplify(get(diagrams, diagram, zero(D)) + value)
            if iszero(combined)
                delete!(diagrams, diagram)
            else
                diagrams[diagram] = combined
            end
        end
    end

    return TwoPIEffectiveAction{D,S,O,E,E2}(diagrams, parameters(L)^O, field_families(L))
end
