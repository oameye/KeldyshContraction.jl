# GraphCombinations backend for symmetry-aware Wick matching.
#
# This layer is intentionally internal. GraphCombinations owns generic colored-port matching,
# automorphism quotienting, and exact integer multiplicities. KeldyshContraction retains all
# field semantics, local contraction admissibility, physical Wick filters, simplification signs,
# and final diagram canonicalization.

@inline _gc_normalized_field(field::Field{S}) where {S<:Statistics} =
    reconstruct(field; position=Bulk(1))

function _gc_intern_field_color!(
    colors::Vector{Field{S}}, field::Field{S}
) where {S<:Statistics}
    normalized = _gc_normalized_field(field)
    index = findfirst(x -> isequal(x, normalized), colors)
    if isnothing(index)
        push!(colors, normalized)
        return length(colors)
    end
    return index
end

function _gc_matching_positions(args_nc::Vector{Field{S}}) where {S<:Statistics}
    positions = Position[]
    for field in args_nc
        p = position(field)
        any(existing -> isequal(existing, p), positions) || push!(positions, p)
    end
    sort!(positions; by=p -> if is_in(p)
        (0, 0)
    elseif is_out(p)
        (1, 0)
    else
        (2, Int(p))
    end)
    return positions
end

function _gc_vertex_index(
    positions::Vector{Position}, field::Field{S}
) where {S<:Statistics}
    index = findfirst(p -> isequal(p, position(field)), positions)
    isnothing(index) &&
        error("Internal error: Wick field position is missing from port vertices.")
    return index
end

function _gc_vertex_signature(args_nc::Vector{Field{S}}, p::Position) where {S<:Statistics}
    signature = Field{S}[
        _gc_normalized_field(field) for field in args_nc if isequal(position(field), p)
    ]
    sort!(signature)
    return signature
end

function _gc_vertex_colors(
    args_nc::Vector{Field{S}}, positions::Vector{Position}
) where {S<:Statistics}
    num_fixed = count(p -> is_in(p) || is_out(p), positions)
    colors = zeros(Int, length(positions))
    signatures = Vector{Vector{Field{S}}}()
    for vertex in eachindex(positions)
        if vertex <= num_fixed
            colors[vertex] = vertex
            continue
        end
        signature = _gc_vertex_signature(args_nc, positions[vertex])
        index = findfirst(existing -> isequal(existing, signature), signatures)
        if isnothing(index)
            push!(signatures, signature)
            index = length(signatures)
        end
        colors[vertex] = num_fixed + index
    end
    return colors, num_fixed
end

function _gc_build_wick_problem(
    args_nc::Vector{Field{S}}, ::Val{E}; regularise=true, _set_reg_to_zero=false
) where {S<:Statistics,E}
    destroys, creates = prepare_args(args_nc, Val(E))
    positions = _gc_matching_positions(args_nc)
    colors, num_fixed = _gc_vertex_colors(args_nc, positions)

    source_colors = Field{S}[]
    target_colors = Field{S}[]
    source_color = Vector{Int}(undef, E)
    target_color = Vector{Int}(undef, E)
    source_vertex = Vector{Int}(undef, E)
    target_vertex = Vector{Int}(undef, E)

    @inbounds for k in 1:E
        source_color[k] = _gc_intern_field_color!(source_colors, destroys[k])
        target_color[k] = _gc_intern_field_color!(target_colors, creates[k])
        source_vertex[k] = _gc_vertex_index(positions, destroys[k])
        target_vertex[k] = _gc_vertex_index(positions, creates[k])
    end

    source_ports = zeros(Int, length(positions), length(source_colors))
    target_ports = zeros(Int, length(positions), length(target_colors))
    @inbounds for k in 1:E
        source_ports[source_vertex[k], source_color[k]] += 1
        target_ports[target_vertex[k], target_color[k]] += 1
    end

    ps = map(position, args_nc)
    candidates = wick_candidates(
        destroys,
        creates,
        Val(E);
        regularise,
        _set_reg_to_zero,
        skip_external_pair=has_in(ps) && has_out(ps),
    )

    compatibility = falses(
        length(positions), length(source_colors), length(positions), length(target_colors)
    )
    lookup = Dict{NTuple{4,Int},Contraction{S}}()
    @inbounds for k in 1:E
        source_cell = (source_vertex[k], source_color[k])
        for (l, contraction) in candidates[k]
            cell = (source_cell[1], source_cell[2], target_vertex[l], target_color[l])
            compatibility[cell...] = true
            if haskey(lookup, cell)
                isequal(lookup[cell], contraction) || error(
                    "Internal error: colored-port encoding merged distinct contractions.",
                )
            else
                lookup[cell] = contraction
            end
        end
    end

    problem = GraphCombinations.ColoredPortProblem(
        colors, source_ports, target_ports, compatibility, num_fixed
    )
    return problem, lookup
end

function _gc_completion_contractions(
    completion::GraphCombinations.WeightedPortCompletion,
    lookup::Dict{NTuple{4,Int},Contraction{S}},
    ::Val{E},
) where {S<:Statistics,E}
    length(completion.edges) == E ||
        error("Internal error: GraphCombinations returned an incomplete Wick matching.")
    contractions = Vector{Contraction{S}}(undef, E)
    @inbounds for i in eachindex(completion.edges)
        edge = completion.edges[i]
        contractions[i] = lookup[(
            edge.source, edge.source_color, edge.target, edge.target_color
        )]
    end
    return contractions
end

function _gc_bosonic_wick_contraction_with_stats(
    args_nc::Vector{Field{Boson}},
    ::Val{E},
    ::Val{E2};
    regularise=true,
    _set_reg_to_zero=false,
    simplify=false,
) where {E,E2}
    problem, lookup = _gc_build_wick_problem(args_nc, Val(E); regularise, _set_reg_to_zero)
    completions, stats = GraphCombinations.generate_weighted_with_stats(problem)

    canonical_weights = Dict{FixedVector{E,Contraction{Boson}},BigInt}()
    for completion in completions
        contractions = _gc_completion_contractions(completion, lookup, Val(E))
        passes_wick_filters(contractions) || continue

        final_contractions, simplification_sign = if simplify
            advanced_to_retarded(contractions, 1)
        else
            contractions, 1
        end
        canonical = canonicalize(final_contractions)
        canonical_key = sorted_wick_key(canonical, Val(E))
        final_weight = completion.weight * Int(simplification_sign)
        canonical_weights[canonical_key] =
            get(canonical_weights, canonical_key, big(0)) + final_weight
    end

    wick_pairings = Tuple{WickPairing{Boson,E},FixedVector{E2,Int},Int}[]
    sizehint!(wick_pairings, length(canonical_weights))
    for (canonical, weight) in canonical_weights
        iszero(weight) && continue
        contractions = Contraction{Boson}[contraction for contraction in canonical]
        topology = legacy_topology(contractions, Val(E2))
        pairing = WickPairing(contractions, Int8(sign(weight)), Val(E))
        multiplicity = Int(abs(weight))
        push!(wick_pairings, (pairing, topology, multiplicity))
    end
    return wick_pairings, stats
end

function _gc_bosonic_wick_contraction(
    args_nc::Vector{Field{Boson}}, ::Val{E}, ::Val{E2}; kwargs...
) where {E,E2}
    pairings, _ = _gc_bosonic_wick_contraction_with_stats(
        args_nc, Val(E), Val(E2); kwargs...
    )
    return pairings
end
