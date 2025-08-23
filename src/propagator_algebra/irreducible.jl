"""
    is_irreducible(multiplicity_list::Vector{Int}, num_vertices::Int) -> Bool

Check if a graph represented by edge multiplicities is irreducible.
A graph is irreducible if removing any single edge does not disconnect it.

# Arguments
- `multiplicity_list`: Vector of integers representing the multiplicity of edges
  between pairs of vertices (in some canonical ordering)
- `num_vertices`: Number of vertices in the graph

# Returns
- `true` if the graph is irreducible (cannot be disconnected by removing one edge)
- `false` if the graph is reducible (can be disconnected by removing one edge)

# Algorithm
1. Check if the original graph is connected using BFS on the multiplicity list
2. Find all edges with multiplicity 1 (potential bridges)
3. For each such edge, temporarily remove it (set multiplicity to 0)
4. Check if the resulting graph is still connected using BFS
5. If any removal disconnects the graph, return false (reducible)
6. If all removals keep the graph connected, return true (irreducible)
"""
function is_irreducible(multiplicity_list::Vector{Int}, num_vertices::Int)
    # If graph has fewer than 2 vertices, it's trivially irreducible
    if num_vertices < 2
        return true
    end

    # If the original graph is not connected, it's not irreducible
    if !is_connected(multiplicity_list, num_vertices)
        return false
    end

    # Find all positions with multiplicity 1 (potential bridges)
    bridge_candidates = findall(isone, multiplicity_list)

    # Test removing each edge with multiplicity 1
    for bridge_idx in bridge_candidates
        # Create a copy of multiplicity list with this edge removed
        test_multiplicity = copy(multiplicity_list)
        test_multiplicity[bridge_idx] = 0

        # If removing this edge disconnects the graph, it's reducible
        if !is_connected(test_multiplicity, num_vertices)
            return false
        end
    end

    # If we reach here, removing any single edge keeps the graph connected
    return true
end

# Check connectivity directly from multiplicity list using BFS
function is_connected(mult_list::Vector{Int}, n_vertices::Int)
    if n_vertices <= 1
        return true
    end

    # Helper function to get multiplicity between vertices i and j
    function get_multiplicity(i::Int, j::Int)
        if i == j
            return 0  # No self-loops
        end
        if i > j
            i, j = j, i  # Ensure i < j for canonical ordering
        end

        # Calculate index in multiplicity list for edge (i,j)
        # For vertices numbered 1 to n, edge (i,j) with i<j maps to:
        # Index = (i-1)*(2n-i)/2 + (j-i)
        # This is the standard upper triangular matrix indexing
        idx = (i - 1) * (2 * n_vertices - i) ÷ 2 + (j - i)

        return (idx <= length(mult_list)) ? mult_list[idx] : 0
    end

    # BFS to check connectivity
    visited = Set{Int}()
    queue = [1]  # Start from vertex 1
    push!(visited, 1)

    while !isempty(queue)
        current = popfirst!(queue)

        # Check all other vertices for connections
        for neighbor in 1:n_vertices
            if neighbor != current && neighbor ∉ visited
                # Check if there's an edge between current and neighbor
                if get_multiplicity(current, neighbor) > 0
                    push!(queue, neighbor)
                    push!(visited, neighbor)
                end
            end
        end
    end

    # Graph is connected if we visited all vertices
    return length(visited) == n_vertices
end
