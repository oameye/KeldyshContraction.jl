using KeldyshContraction
using Combinatorics
using GraphCombinations
using SmallCollections: FixedVector

import GraphCombinations as GC
import KeldyshContraction as KC

struct AdapterStats
    terms::Int
    automorphisms::Int
    quotient_states::Int
    transitions::Int
    merged_transitions::Int
    completed_states::Int
    filter_survivors::Int
    canonical_unique::Int
end
AdapterStats() = AdapterStats(0, 0, 0, 0, 0, 0, 0, 0)
Base.:+(a::AdapterStats, b::AdapterStats) = AdapterStats(
    a.terms + b.terms,
    a.automorphisms + b.automorphisms,
    a.quotient_states + b.quotient_states,
    a.transitions + b.transitions,
    a.merged_transitions + b.merged_transitions,
    a.completed_states + b.completed_states,
    a.filter_survivors + b.filter_survivors,
    a.canonical_unique + b.canonical_unique,
)

normalized_field(field::KC.Field) = KC.reconstruct(field; position=KC.Bulk(1))

function intern_field_color!(colors::Vector{KC.Field{S}}, field::KC.Field{S}) where {S}
    normalized = normalized_field(field)
    index = findfirst(x -> isequal(x, normalized), colors)
    if isnothing(index)
        push!(colors, normalized)
        return length(colors)
    end
    return index
end

function matching_positions(args_nc)
    positions = KC.Position[]
    for field in args_nc
        position = KC.position(field)
        any(p -> isequal(p, position), positions) || push!(positions, position)
    end
    sort!(positions; by=p -> KC.is_in(p) ? (0, 0) : KC.is_out(p) ? (1, 0) : (2, Int(p)))
    return positions
end

function vertex_index(positions, field)
    index = findfirst(p -> isequal(p, KC.position(field)), positions)
    isnothing(index) && error("missing field position")
    return index
end

function vertex_signature(args_nc::Vector{KC.Field{S}}, position) where {S}
    signature = KC.Field{S}[
        normalized_field(field) for field in args_nc if isequal(KC.position(field), position)
    ]
    sort!(signature)
    return signature
end

function vertex_colors(args_nc::Vector{KC.Field{S}}, positions) where {S}
    num_fixed = count(p -> KC.is_in(p) || KC.is_out(p), positions)
    colors = zeros(Int, length(positions))
    signatures = Vector{Vector{KC.Field{S}}}()
    for vertex in eachindex(positions)
        if vertex <= num_fixed
            colors[vertex] = vertex
            continue
        end
        signature = vertex_signature(args_nc, positions[vertex])
        index = findfirst(x -> isequal(x, signature), signatures)
        if isnothing(index)
            push!(signatures, signature)
            index = length(signatures)
        end
        colors[vertex] = num_fixed + index
    end
    return colors, num_fixed
end

function build_gc_problem(
    args_nc::Vector{KC.Field{S}}, ::Val{E}; regularise=true, _set_reg_to_zero=false
) where {S<:KC.Statistics,E}
    destroys, creates = KC.prepare_args(args_nc, Val(E))
    positions = matching_positions(args_nc)
    colors, num_fixed = vertex_colors(args_nc, positions)
    source_colors, target_colors = KC.Field{S}[], KC.Field{S}[]
    source_color, target_color = Vector{Int}(undef, E), Vector{Int}(undef, E)
    source_vertex, target_vertex = Vector{Int}(undef, E), Vector{Int}(undef, E)
    for k in 1:E
        source_color[k] = intern_field_color!(source_colors, destroys[k])
        target_color[k] = intern_field_color!(target_colors, creates[k])
        source_vertex[k] = vertex_index(positions, destroys[k])
        target_vertex[k] = vertex_index(positions, creates[k])
    end

    source_ports = zeros(Int, length(positions), length(source_colors))
    target_ports = zeros(Int, length(positions), length(target_colors))
    for k in 1:E
        source_ports[source_vertex[k], source_color[k]] += 1
        target_ports[target_vertex[k], target_color[k]] += 1
    end

    ps = map(KC.position, args_nc)
    candidates = KC.wick_candidates(
        destroys,
        creates,
        Val(E);
        regularise,
        _set_reg_to_zero,
        skip_external_pair=KC.has_in(ps) && KC.has_out(ps),
    )
    compatibility = falses(
        length(positions), length(source_colors), length(positions), length(target_colors)
    )
    lookup = Dict{NTuple{4,Int},KC.Contraction{S}}()
    for k in 1:E
        sv, sc = source_vertex[k], source_color[k]
        for (l, contraction) in candidates[k]
            cell = (sv, sc, target_vertex[l], target_color[l])
            compatibility[cell...] = true
            if haskey(lookup, cell)
                isequal(lookup[cell], contraction) ||
                    error("port coloring merged distinct KC contractions")
            else
                lookup[cell] = contraction
            end
        end
    end
    problem = GC._PortMatchingProblem(
        colors, source_ports, target_ports, compatibility, num_fixed
    )
    return problem, lookup
end

function expected_weights(
    args_nc::Vector{KC.Field{S}},
    ::Val{E},
    ::Val{E2};
    regularise=true,
    _set_reg_to_zero=false,
    simplify=false,
) where {S<:KC.Statistics,E,E2}
    pairings = KC._wick_contraction(
        args_nc, Val(E), Val(E2); regularise, _set_reg_to_zero, simplify
    )
    weights = Dict{FixedVector{E,KC.Contraction{S}},BigInt}()
    for (pairing, _, multiplicity) in pairings
        weights[pairing.contractions] =
            get(weights, pairing.contractions, big(0)) +
            big(Int(pairing.sign)) * multiplicity
    end
    filter!(pair -> !iszero(last(pair)), weights)
    return weights
end

function adapter_weights(
    args_nc::Vector{KC.Field{S}},
    ::Val{E};
    regularise=true,
    _set_reg_to_zero=false,
    simplify=false,
) where {S<:KC.Statistics,E}
    problem, lookup = build_gc_problem(args_nc, Val(E); regularise, _set_reg_to_zero)
    results, search = GC._weighted_port_matchings_with_stats(problem)
    weights = Dict{FixedVector{E,KC.Contraction{S}},BigInt}()
    filter_survivors = 0
    for (edges, multiplicity) in results
        contractions = Vector{KC.Contraction{S}}(undef, E)
        for i in eachindex(edges)
            edge = edges[i]
            contractions[i] = lookup[
                (edge.source, edge.source_color, edge.target, edge.target_color)
            ]
        end
        KC.passes_wick_filters(contractions) || continue
        filter_survivors += 1
        final_contractions, sign = simplify ?
            KC.advanced_to_retarded(contractions, 1) : (contractions, 1)
        key = KC.sorted_wick_key(KC.canonicalize(final_contractions), Val(E))
        weights[key] = get(weights, key, big(0)) + multiplicity * Int(sign)
    end
    filter!(pair -> !iszero(last(pair)), weights)
    stats = AdapterStats(
        1,
        search.automorphisms,
        sum(search.layer_states),
        search.transitions,
        search.merged_transitions,
        length(results),
        filter_survivors,
        length(weights),
    )
    return weights, stats
end

function validate_fields(
    args_nc::Vector{KC.Field{Boson}}, ::Val{E}, ::Val{E2}; kwargs...
) where {E,E2}
    expected = expected_weights(args_nc, Val(E), Val(E2); kwargs...)
    actual, stats = adapter_weights(args_nc, Val(E); kwargs...)
    actual == expected || error("GC adapter disagrees with the KC bosonic Wick oracle")
    return stats
end

function validate_component(
    in_out::KC.QMul{CI,Boson},
    L::KC.InteractionLagrangian{CL,Boson},
    ::Val{O},
    ::Val{E};
    simplify=true,
    _set_reg_to_zero=true,
) where {CI<:Number,CL<:Number,O,E}
    total = AdapterStats()
    regularise = KC.should_regularise(L.lagrangian)
    for coefficients in Combinatorics.multiexponents(length(L.lagrangian), O)
        idxs = KC.indices_from_counts(coefficients)
        mult = Combinatorics.multinomial(coefficients...)
        qmul = mult * prod(L(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))
        total += validate_fields(
            (in_out * qmul).args_nc,
            Val(E),
            Val(KC.max_edges(O));
            regularise,
            _set_reg_to_zero,
            simplify,
        )
    end
    return total
end

function print_stats(workload, component, s::AdapterStats)
    println(
        join(
            (
                workload,
                component,
                s.terms,
                s.automorphisms,
                s.quotient_states,
                s.transitions,
                s.merged_transitions,
                s.completed_states,
                s.filter_survivors,
                s.canonical_unique,
            ),
            '\t',
        ),
    )
end

function validate_interaction(
    workload,
    L::KC.InteractionLagrangian{C,Boson},
    order::Int,
    edges::Int;
    simplify=true,
    _set_reg_to_zero=true,
) where {C<:Number}
    products = KC.propagator_external_products(Boson, KC.propagator_fields(L, nothing)...)
    total = AdapterStats()
    for (name, in_out) in zip(("K", "R", "A"), products)
        stats = validate_component(
            in_out, L, Val(order), Val(edges); simplify, _set_reg_to_zero
        )
        print_stats(workload, name, stats)
        total += stats
    end
    print_stats(workload, "total", total)
    return total
end

println(
    "workload\tcomponent\tterms\tautomorphisms\tquotient_states\ttransitions\t",
    "merged_transitions\tcompleted_states\tfilter_survivors\tcanonical_unique",
)

@qfields ϕadapter::Boson
c, q = ϕadapter[Classical], ϕadapter[Quantum]
elastic = -(
    0.5 * (c^2 + q^2) * bar(c) * bar(q) +
    0.5 * c * q * (bar(c)^2 + bar(q)^2)
)
L_g = InteractionLagrangian(elastic, :g)
loss =
    0.5 * bar(c) * bar(q) *
    (c(KC.Regularisation.Minus)^2 + q(KC.Regularisation.Minus)^2) -
    0.5 * c(KC.Regularisation.Plus) * q(KC.Regularisation.Plus) *
    (bar(c)^2 + bar(q)^2) +
    bar(c) * bar(q) *
    (c(KC.Regularisation.Plus) * q(KC.Regularisation.Plus) +
     c(KC.Regularisation.Minus) * q(KC.Regularisation.Minus))
L_γ = InteractionLagrangian(loss, :γ)

validate_interaction("boson_g2", L_g, 2, 5)
validate_interaction("boson_gamma2", L_γ, 2, 5; simplify=true, _set_reg_to_zero=true)
validate_interaction("boson_g3", L_g, 3, 7)
