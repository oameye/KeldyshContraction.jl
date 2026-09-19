# GraphCombinations backend for symmetry-aware Wick generation.
#
# GC owns generic colored-port residual-count generation, automorphism quotienting,
# canonical relabeling, and exact integer accumulation. KC retains field semantics,
# fermionic orientation/sign transport, physical filters, final propagator
# canonicalization, and historical topology.

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
    result = Position[]
    for field in args_nc
        p = position(field)
        any(existing -> isequal(existing, p), result) || push!(result, p)
    end
    sort!(result; by=p -> if is_out(p)
        (0, 0)
    elseif is_in(p)
        (1, 0)
    else
        (2, Int(p))
    end)
    return result
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
    num_fixed = count(p -> is_out(p) || is_in(p), positions)
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

struct _GCFermionPortTransport
    source_vertices::Vector{Int}
    source_colors::Vector{Int}
    target_vertices::Vector{Int}
    target_colors::Vector{Int}
end

@inline function _gc_cell_isless(a::NTuple{2,Int}, b::NTuple{2,Int})
    a[1] == b[1] || return a[1] < b[1]
    return a[2] < b[2]
end

function _gc_order_sign(cells::Vector{NTuple{2,Int}})::Int
    sign = 1
    @inbounds for i in 1:(length(cells) - 1)
        for j in (i + 1):length(cells)
            _gc_cell_isless(cells[j], cells[i]) && (sign = -sign)
        end
    end
    return sign
end

function _gc_mapped_order_sign(
    cells::Vector{NTuple{2,Int}}, witness::GC.PortRelabeling
)::Int
    mapped = Vector{NTuple{2,Int}}(undef, length(cells))
    @inbounds for i in eachindex(cells)
        vertex, color = cells[i]
        mapped[i] = (witness.vertex_map[vertex], color)
    end
    return _gc_order_sign(mapped)
end

function _gc_expanded_cells(
    counts::AbstractMatrix{<:Integer}, added_vertex::Int=0, added_color::Int=0
)
    extra = Int(!iszero(added_vertex))
    cells = Vector{NTuple{2,Int}}()
    sizehint!(cells, sum(counts) + extra)
    @inbounds for vertex in axes(counts, 1)
        for color in axes(counts, 2)
            multiplicity = Int(counts[vertex, color])
            if vertex == added_vertex && color == added_color
                multiplicity += 1
            end
            for _ in 1:multiplicity
                push!(cells, (vertex, color))
            end
        end
    end
    return cells
end

function GC.initial_port_weight(
    transport::_GCFermionPortTransport,
    ::GC.ColoredPortState,
    ::GC.ColoredPortState,
    witness::GC.PortRelabeling,
)::BigInt
    source_cells = Vector{NTuple{2,Int}}(undef, length(transport.source_vertices))
    target_cells = Vector{NTuple{2,Int}}(undef, length(transport.target_vertices))
    @inbounds for i in eachindex(source_cells)
        source_cells[i] = (
            witness.vertex_map[transport.source_vertices[i]], transport.source_colors[i]
        )
        target_cells[i] = (
            witness.vertex_map[transport.target_vertices[i]], transport.target_colors[i]
        )
    end
    return big(_gc_order_sign(source_cells) * _gc_order_sign(target_cells))
end

function GC.transport_port_weight(
    ::_GCFermionPortTransport,
    parent_weight::BigInt,
    multiplicity::Int,
    added_edge::GC.ColoredPortEdge,
    raw_child::GC.ColoredPortState,
    ::GC.ColoredPortState,
    witness::GC.PortRelabeling,
)::BigInt
    source_cell = (added_edge.source, added_edge.source_color)
    target_cell = (added_edge.target, added_edge.target_color)

    parent_sources = _gc_expanded_cells(
        GC.source_port_counts(raw_child), added_edge.source, added_edge.source_color
    )
    parent_targets = _gc_expanded_cells(
        GC.target_port_counts(raw_child), added_edge.target, added_edge.target_color
    )

    source_rank = findfirst(==(source_cell), parent_sources)
    isnothing(source_rank) && error("Internal error: selected Wick source slot is missing.")
    target_ranks = findall(==(target_cell), parent_targets)
    length(target_ranks) == multiplicity || error(
        "Internal error: GC target multiplicity disagrees with the concrete fermionic slot count.",
    )

    residual_sources = copy(parent_sources)
    deleteat!(residual_sources, source_rank)
    source_relabel_sign = _gc_mapped_order_sign(residual_sources, witness)

    signed_minor = 0
    for target_rank in target_ranks
        residual_targets = copy(parent_targets)
        deleteat!(residual_targets, target_rank)
        target_relabel_sign = _gc_mapped_order_sign(residual_targets, witness)
        cofactor_sign = isodd(source_rank + target_rank) ? -1 : 1
        signed_minor += cofactor_sign * source_relabel_sign * target_relabel_sign
    end
    return parent_weight * signed_minor
end

function _gc_build_wick_problem(
    args_nc::Vector{Field{S}}, ::Val{E}; regularise=true, _set_reg_to_zero=false
) where {S<:Statistics,E}
    destroys, creates = prepare_args(args_nc, Val(E))
    positions = _gc_matching_positions(args_nc)
    colors, num_fixed = _gc_vertex_colors(args_nc, positions)

    source_fields = Field{S}[]
    target_fields = Field{S}[]
    source_color = Vector{Int}(undef, E)
    target_color = Vector{Int}(undef, E)
    source_vertex = Vector{Int}(undef, E)
    target_vertex = Vector{Int}(undef, E)

    @inbounds for k in 1:E
        source_color[k] = _gc_intern_field_color!(source_fields, destroys[k])
        target_color[k] = _gc_intern_field_color!(target_fields, creates[k])
        source_vertex[k] = _gc_vertex_index(positions, destroys[k])
        target_vertex[k] = _gc_vertex_index(positions, creates[k])
    end

    source_ports = zeros(Int, length(positions), length(source_fields))
    target_ports = zeros(Int, length(positions), length(target_fields))
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
        length(positions), length(source_fields), length(positions), length(target_fields)
    )
    lookup = Dict{NTuple{4,Int},Contraction{S}}()
    @inbounds for k in 1:E
        for (l, contraction) in candidates[k]
            cell = (source_vertex[k], source_color[k], target_vertex[l], target_color[l])
            compatibility[cell...] = true
            if haskey(lookup, cell)
                isequal(lookup[cell], contraction) || error(
                    "Internal error: colored-port encoding merged distinct Wick contractions.",
                )
            else
                lookup[cell] = contraction
            end
        end
    end

    problem = GC.ColoredPortProblem(
        colors, source_ports, target_ports, compatibility, num_fixed
    )
    transport = if S <: Fermion
        _GCFermionPortTransport(source_vertex, source_color, target_vertex, target_color)
    else
        GC.MultiplicityPortTransport()
    end
    return problem, lookup, transport
end

function _gc_completion_contractions(
    completion::GC.WeightedPortCompletion,
    lookup::Dict{NTuple{4,Int},Contraction{S}},
    ::Val{E},
) where {S<:Statistics,E}
    length(completion.edges) == E ||
        error("Internal error: GraphCombinations returned an incomplete Wick matching.")
    contractions = Vector{Contraction{S}}(undef, E)
    @inbounds for i in eachindex(completion.edges)
        edge = completion.edges[i]
        key = (edge.source, edge.source_color, edge.target, edge.target_color)
        haskey(lookup, key) || error(
            "Internal error: canonical GC Wick edge has no physical contraction representative.",
        )
        contractions[i] = lookup[key]
    end
    return contractions
end

function _gc_wick_contraction_with_stats(
    args_nc::Vector{Field{S}},
    ::Val{E},
    ::Val{E2};
    regularise=true,
    _set_reg_to_zero=false,
    simplify=false,
) where {S<:Statistics,E,E2}
    problem, lookup, transport = _gc_build_wick_problem(
        args_nc, Val(E); regularise, _set_reg_to_zero
    )
    completions, stats = GC.generate_weighted_with_stats(problem; transport)

    canonical_scratch = physical_canonicalization_workspace(Val(E))
    canonical_weights = Dict{FixedVector{E,Contraction{S}},BigInt}()
    for completion in completions
        contractions = _gc_completion_contractions(completion, lookup, Val(E))
        passes_wick_filters(contractions) || continue

        final_contractions, simplification_sign = if simplify
            advanced_to_retarded(contractions, 1)
        else
            contractions, 1
        end
        canonical = canonicalize(final_contractions, canonical_scratch)
        canonical_key = sorted_wick_key(canonical, Val(E))
        final_weight = completion.weight * Int(simplification_sign)
        canonical_weights[canonical_key] =
            get(canonical_weights, canonical_key, big(0)) + final_weight
    end

    wick_pairings = Tuple{WickPairing{S,E},FixedVector{E2,Int},Int}[]
    sizehint!(wick_pairings, length(canonical_weights))
    for (canonical, weight) in canonical_weights
        iszero(weight) && continue
        contractions = Contraction{S}[contraction for contraction in canonical]
        topology = legacy_topology(contractions, Val(E2))
        pairing = WickPairing(contractions, Int8(sign(weight)), Val(E))
        push!(wick_pairings, (pairing, topology, Int(abs(weight))))
    end
    return wick_pairings, stats
end

function _gc_wick_contraction(
    args_nc::Vector{Field{S}}, ::Val{E}, ::Val{E2}; kwargs...
) where {S<:Statistics,E,E2}
    pairings, _ = _gc_wick_contraction_with_stats(args_nc, Val(E), Val(E2); kwargs...)
    return pairings
end
