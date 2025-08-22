"""
    advanced_to_retarded(x::T) where {T<:SymbolicUtils.Symbolic}

Apply the transformation to change the advanced propagator to retarded:

``G^A(y, y)=-G^R(y, y)``

with ``y =(\\vec{y},t)``.
Note the expression is only valid for equal space-time coordinates.
"""
function advanced_to_retarded(
    x::Vector{Contraction}, prefactor::Number
)::Tuple{Vector{Contraction},Number}
    ff(x::Contraction) = is_advanced(x) && same_position(x)
    adv_idx = findall(ff, x)
    if isempty(adv_idx)
        return x, prefactor
    end
    x′ = deepcopy(x)
    # TODO or make x immutable or change in place
    for i in adv_idx
        prefactor *= -1
        x′[i] = adjoint(x[i])
    end
    return x′, prefactor
end

function sort_by_position_and_type(p::Contraction)::Float64
    if is_out(p)
        return -Inf
    elseif is_in(p)
        return Inf
    else
        i, j = integer_positions(p)
        type = Int(propagator_type(p...))
        return float(pairing(i, j) * 3 + type)
    end
end
sort_by_position_and_type(p::Edge)::Float64 = sort_by_position_and_type(fields(p))

function canonicalize(vs::Vector{Contraction})
    canonical = is_canonical(vs)
    if !canonical
        # Extract positions and canonicalize
        position_edges = [positions(c) for c in vs]
        canonical_position_edges = canonicalize(position_edges)

        # Apply the relabeling back to the contractions
        vs = map(zip(vs, canonical_position_edges)) do (c, new_positions)
            # c is a tuple of fields (Destroy, Create)
            # new_positions is a tuple of positions (Position, Position)
            field1, field2 = c
            pos1, pos2 = new_positions
            (field1(pos1), field2(pos2))
        end
    end
    return vs
end

function canonicalize(position_edges::Vector{NTuple{2, Position}})
    canonical = is_canonical(position_edges)
    if !canonical
        bulk_mapping = find_canonical_bulk_mapping(position_edges)
        position_edges = apply_bulk_relabeling(position_edges, bulk_mapping)
    end

    # Sort edges in canonical order for complete canonical form
    # Sort by: (position1_index, position2_index) where external positions have special indices
    sorted_edges = sort(position_edges, by = edge -> (Int(edge[1]), Int(edge[2])))

    return sorted_edges
end

"""
Create a canonical ordering signature for a node based on its structural properties.
This includes recursive neighborhood information for better discrimination.
"""
function node_signature(node::Int, adjacency::Dict{Int, Set{Int}}, all_bulk_positions::Set{Int}, depth::Int=2)::Vector{Int}
    if depth == 0
        return [0]
    end

    neighbors = get(adjacency, node, Set{Int}())
    bulk_neighbors = [n for n in neighbors if n ∈ all_bulk_positions]

    # Start with degree
    signature = [length(bulk_neighbors)]

    # Add sorted neighbor signatures recursively
    neighbor_sigs = []
    for neighbor in sort(bulk_neighbors)
        if depth > 1
            neighbor_sig = node_signature(neighbor, adjacency, all_bulk_positions, depth-1)
            push!(neighbor_sigs, neighbor_sig)
        else
            push!(neighbor_sigs, [length(get(adjacency, neighbor, Set{Int}()))])
        end
    end

    # Sort neighbor signatures and flatten
    sort!(neighbor_sigs)
    for sig in neighbor_sigs
        append!(signature, sig)
    end

    return signature
end

"""
Find canonical bulk node mapping using full graph isomorphism for unique labeling.
This version considers the complete graph structure including external connections.
"""
function find_canonical_bulk_mapping(position_edges::Vector{NTuple{2, Position}})::Dict{Int, Int}
    # Get all bulk positions in the graph
    all_bulk_positions = get_all_bulk_positions(position_edges)

    if length(all_bulk_positions) <= 1
        # For single or no bulk positions, map to canonical form
        canonical_mapping = Dict{Int, Int}()
        if length(all_bulk_positions) == 1
            old_idx = first(all_bulk_positions)
            canonical_mapping[old_idx] = 1
        end
        return canonical_mapping
    end

    # Build adjacency list from position edges (including external connections)
    adjacency = build_adjacency_list(position_edges)

    # Find the bulk node connected to Out (this will always be labeled as 1)
    output_connected_bulk = Set{Int}()
    for edge in position_edges
        if edge[1] == Out() && is_bulk(edge[2])
            push!(output_connected_bulk, index(edge[2]))
        elseif edge[2] == Out() && is_bulk(edge[1])
            push!(output_connected_bulk, index(edge[1]))
        end
    end

    # Find the bulk node connected to In (this helps break symmetry)
    input_connected_bulk = Set{Int}()
    for edge in position_edges
        if edge[1] == In() && is_bulk(edge[2])
            push!(input_connected_bulk, index(edge[2]))
        elseif edge[2] == In() && is_bulk(edge[1])
            push!(input_connected_bulk, index(edge[1]))
        end
    end

    # Start canonical labeling
    canonical_mapping = Dict{Int, Int}()
    visited = Set{Int}()
    next_canonical_label = 1

    # Always label the output-connected node as 1
    if length(output_connected_bulk) == 1
        output_node = first(output_connected_bulk)
        canonical_mapping[output_node] = next_canonical_label
        next_canonical_label += 1
        push!(visited, output_node)
    end

    # For remaining nodes, use a canonical ordering based on:
    # 1. Distance from output node (BFS layers)
    # 2. Structural signature within each layer
    # 3. Special handling for input-connected nodes

    current_level = [node for node in output_connected_bulk if node ∉ visited]

    while !isempty(current_level) || length(visited) < length(all_bulk_positions)
        if isempty(current_level)
            # Handle disconnected components
            current_level = [node for node in all_bulk_positions if node ∉ visited]
        end

        # Create enhanced signatures that include external connectivity
        node_signatures = []
        for node in current_level
            if node ∉ visited
                signature = node_signature(node, adjacency, all_bulk_positions)

                # Add external connectivity information to signature
                is_input_connected = node ∈ input_connected_bulk
                is_output_connected = node ∈ output_connected_bulk

                # Enhanced signature: (base_signature, input_connected, output_connected, node_index)
                enhanced_signature = (signature, is_input_connected, is_output_connected, node)
                push!(node_signatures, (node, enhanced_signature))
            end
        end

        # Sort by enhanced signature for canonical ordering
        sort!(node_signatures, by=x -> x[2])

        next_level = Int[]

        for (node, _) in node_signatures
            if node ∉ visited
                canonical_mapping[node] = next_canonical_label
                next_canonical_label += 1
                push!(visited, node)

                # Add unvisited neighbors to next level
                if haskey(adjacency, node)
                    for neighbor in adjacency[node]
                        if neighbor ∈ all_bulk_positions && neighbor ∉ visited && neighbor ∉ next_level
                            push!(next_level, neighbor)
                        end
                    end
                end
            end
        end

        current_level = next_level
    end

    return canonical_mapping
end

"""
Build adjacency list representation of the graph from position edges.
"""
function build_adjacency_list(position_edges::Vector{NTuple{2, Position}})::Dict{Int, Set{Int}}
    adjacency = Dict{Int, Set{Int}}()

    for edge in position_edges
        pos1, pos2 = edge
        idx1, idx2 = index(pos1), index(pos2)

        # Add edges (skip external nodes for adjacency purposes)
        if is_bulk(pos1)
            if !haskey(adjacency, idx1)
                adjacency[idx1] = Set{Int}()
            end
            if is_bulk(pos2)
                push!(adjacency[idx1], idx2)
            end
        end

        if is_bulk(pos2)
            if !haskey(adjacency, idx2)
                adjacency[idx2] = Set{Int}()
            end
            if is_bulk(pos1)
                push!(adjacency[idx2], idx1)
            end
        end
    end

    return adjacency
end

"""
Get all bulk position indices present in the position edges.
"""
function get_all_bulk_positions(position_edges::Vector{NTuple{2, Position}})::Set{Int}
    bulk_positions = Set{Int}()
    for edge in position_edges
        for pos in edge
            if is_bulk(pos)
                push!(bulk_positions, index(pos))
            end
        end
    end
    return bulk_positions
end

"""
Apply bulk relabeling to position edges based on the mapping.
"""
function apply_bulk_relabeling(position_edges::Vector{NTuple{2, Position}}, mapping::Dict{Int, Int})
    return map(position_edges) do edge
        map(edge) do pos
            if is_bulk(pos)
                old_idx = index(pos)
                new_idx = get(mapping, old_idx, old_idx)  # Default to old if not in mapping
                Bulk(Int(new_idx))  # Ensure new_idx is Int type
            else
                pos
            end
        end
    end
end

"""
A canonical contraction (diagram) is one where:
1. The output field is connected to bulk position with index 1 (if any bulk positions exist)
2. Bulk positions are labeled in BFS order starting from the output node
3. Bulk positions form a consecutive sequence starting from 1
"""
function is_canonical(vs::Vector{Contraction})
    position_edges = [positions(c) for c in vs]
    return is_canonical(position_edges)
end

function is_canonical(position_edges::Vector{NTuple{2, Position}})
    if isempty(position_edges)
        return true
    end

    # Get all bulk positions
    all_bulk_positions = get_all_bulk_positions(position_edges)

    if isempty(all_bulk_positions)
        return true  # No bulk positions, trivially canonical
    end

    # Check if bulk positions form consecutive sequence starting from 1
    expected_positions = Set(1:length(all_bulk_positions))
    if all_bulk_positions != expected_positions
        return false
    end

    # Check if output is connected to Bulk(1)
    output_edge = position_edges[1]  # First edge contains output
    output_bulk_positions = [pos for pos in output_edge if is_bulk(pos)]

    if Bulk(1) ∉ output_bulk_positions
        return false
    end

    # Check if the current labeling matches the canonical BFS ordering
    canonical_mapping = find_canonical_bulk_mapping(position_edges)

    # If already canonical, mapping should be identity for all positions
    for pos in all_bulk_positions
        if get(canonical_mapping, pos, pos) != pos
            return false
        end
    end

    return true
end
