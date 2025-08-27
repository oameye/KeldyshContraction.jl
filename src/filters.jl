#################################
#       Consistency checks
#################################

is_qq_contraction(v::Contraction) = iszero(sum(Int.(contour.(v))))
# has_qq_contraction(vv::Vector{Vector{<:QField}}) = any(is_qq_contraction.(vv))

function is_physical_propagator(a::Contraction)
    len = length(a) == 2 # propagator has length 2
    positions = position.(a)
    # ∨ Can't make a propagator with In and Out coordinate
    in_out = has_in(positions) ? !has_out(positions) : true
    physical = all(is_physical.(a))
    return len && in_out && physical
end

function contraction_filter(v)
    # istwo = all(length.(v) .== 2) # only two-point contractions
    # if !istwo
    #     return false
    if !is_conserved(v) # contractions with creation/annihilation
        return false
    elseif !is_physical_propagator(v) # propagators are physical
        return false
    else # there should be no quantum-quantum contractions
        return !is_qq_contraction(v)
    end
end

"""
Filter the zero loops from the list of contractions

Gᴿ(1,2) Gᴿ(2,1) = 0
Gᴬ(1,2) Gᴬ(2,1) = 0
Gᴿ(1,2) Gᴬ(1,2) = 0
"""
function has_two_zero_loop(vs::Vector{Contraction})
    ps = positions.(vs)
    sorted_ps = sort_tuple.(ps)
    loops = find_equal_pairs(sorted_ps)
    for loop in loops
        T1 = propagator_type(vs[loop[1]]...)
        T2 = propagator_type(vs[loop[2]]...)
        if all(is_retarded, (T1, T2)) && is_reversed(ps[loop[1]], ps[loop[2]])
            # @info "Has zero loop"
            return true
        elseif all(is_advanced, (T1, T2)) && is_reversed(ps[loop[1]], ps[loop[2]])
            # @info "Has zero loop"
            return true
        elseif isequal(ps[loop[1]], ps[loop[2]]) && retarded_and_advanced_pair(T1, T2)
            # @info "Has zero loop"
            return true
        end
    end
    return false
end
function has_zero_loop(vs::Vector{Contraction})
    if has_two_zero_loop(vs)
        return true
    end

    g, max_label, has_in = make_graph(Graphs.SimpleDiGraph, vs)
    for cycle in Graphs.simplecycles(g)
        if isone(length(cycle))
            continue
        end
        for i in eachindex(cycle)
            cycle[i] -= Int(has_in)
        end

        edges_in_cycle = make_directed_edges(cycle)
        function is_in_cycle(edge)
            ps = integer_positions(edge)
            return ps ∈ edges_in_cycle
        end
        cycle_contractions = filter(is_in_cycle, vs)

        cycle_ptype = map(cc -> propagator_type(cc...), cycle_contractions)
        if count(is_advanced, cycle_ptype) >= length(cycle)
            return true
        elseif count(is_retarded, cycle_ptype) >= length(cycle)
            return true
        end

        # # check for time incostencies $G^A(1,2)G^R(1,3)G^R(3,2)0$
        # # t_1 < t_2 < t_3, t_3 < t_1,
        # cycle_contractions = filter(is_keldysh, cycle_contractions)
        # retarded_idxs = findall(is_retarded, cycle_contractions)
        # for i in retarded_idxs
        #     cycle_contractions[i] = reverse(cycle_contractions[i])
        # end
        # cycle_contractions_ps = map(integer_positions,cycle_contractions)
        # _edges = Graphs.Edge.(cycle_contractions_ps)
        # h_cycle = Graphs.SimpleDiGraph(_edges)
        # if !Graphs.is_connected(h_cycle)
        #     return true
        # end
    end
    return false
end

function make_directed_edges(cycle)
    l = length(cycle)
    return Set((cycle[i], cycle[mod1(i + 1, l)]) for i in eachindex(cycle))
end

function is_reversed(p1, p2)
    return isequal(p1[1], p2[2]) && isequal(p1[2], p2[1])
end
function retarded_and_advanced_pair(T1, T2)
    return (is_retarded(T1) && is_advanced(T2)) || (is_advanced(T1) && is_retarded(T2))
end
"""
    find_equal_pairs(v::Vector)

Find all indices of equal pairs in a vector.

Returns a vector of tuples `(i, j)` where `i < j` and `v[i] == v[j]`.

## Examples
```jldoctest
julia> using KeldyshContraction: find_equal_pairs

julia> find_equal_pairs([1, 2, 1, 3, 2])
2-element Vector{Tuple{Int64, Int64}}:
 (1, 3)
 (2, 5)

julia> find_equal_pairs([1, 2, 3])
Tuple{Int64, Int64}[]

julia> find_equal_pairs([1, 1, 1])
3-element Vector{Tuple{Int64, Int64}}:
 (1, 2)
 (1, 3)
 (2, 3)
```
"""
function find_equal_pairs(vec)
    n = length(vec)
    pairs = Tuple{Int,Int}[]

    for i in 1:n
        for j in (i + 1):n
            if isequal(vec[i], vec[j])
                push!(pairs, (i, j))
            end
        end
    end

    return pairs
end

######################
#     regularise
######################

"""
    regular(p::Average)

Checks if the propagator `p` is regular.
A regular propagator is one that is:
- not at equal time points (tadpole)
- or `p` is not of [`PropagatorType`](@ref) `Retarded` while also having a negative [`Regularisation`](@ref)
- or `p` is not of [`PropagatorType`](@ref) `Advanced` while also having a positive [`Regularisation`](@ref)
"""
function regular(qs::Contraction)
    positions = position.(qs)
    if !isequal(positions[1], positions[2]) # self-contraction
        return true
    end
    _isbulk = is_bulk(qs)
    _reg = regularisations(qs)
    T = propagator_type(qs...)
    if !_isbulk || subtraction(_reg) == 0
        return true
    elseif subtraction(_reg) < 0 && T == PropagatorType.Retarded
        return false
    elseif subtraction(_reg) > 0 && T == PropagatorType.Advanced
        return false
    else
        return true
    end
end

function should_regularise(qmul::QMul)::Bool
    reguralise = false
    for q in qmul.args_nc
        if !iszero(Int(regularisation(q)))
            reguralise = true
            break
        end
    end
    return reguralise
end
function should_regularise(qadd::QAdd)::Bool
    reguralise = false
    for q in qadd.arguments
        if should_regularise(q)
            reguralise = true
            break
        end
    end
    return reguralise
end

function is_connected(vs::Vector{Contraction})
    ps = integer_positions.(vs)
    in_or_out = findfirst(p -> 1 ∈ p || 2 ∈ p, ps) # in case it a vacuum diagram
    edges = isnothing(in_or_out) ? map(p -> p .- 2, ps) : ps
    return is_connected(edges)
end

function is_connected(edges::Union{Vector{Tuple{Int,Int}},Vector{Tuple{Int8,Int8}}})
    # Find all unique vertices
    all_vertices = vertices(edges)
    if isempty(all_vertices)
        return true  # Empty graph is considered connected
    end

    # Find connected components using BFS
    components = connected_components(all_vertices, edges)

    is_single_component = length(components) == 1

    # if !is_single_component
    #     # @info "Contraction is not connected. Found $(length(components)) components:"
    #     for (i, comp) in enumerate(components)
    #         # @info "Component $i: $comp"
    #     end
    # end

    return is_single_component
end

function vertices(ps)
    vertices = Set{Int}()
    for edge in ps
        push!(vertices, edge[1])
        push!(vertices, edge[2])
    end
    return vertices
end

function connected_components(vertices, edges)
    # Find connected components using BFS
    components = Vector{Set{Int}}()
    remaining = Set(vertices)

    while !isempty(remaining)
        component = Set{Int}()
        queue = [first(remaining)]
        push!(component, first(queue))
        delete!(remaining, first(queue))

        while !isempty(queue)
            current = popfirst!(queue)

            # Find neighbors
            for edge in edges
                neighbor = nothing
                if edge[1] == current
                    neighbor = edge[2]
                elseif edge[2] == current
                    neighbor = edge[1]
                end

                if !isnothing(neighbor) && neighbor ∈ remaining
                    push!(queue, neighbor)
                    push!(component, neighbor)
                    delete!(remaining, neighbor)
                end
            end
        end

        push!(components, component)
    end
    return components
end

"""
    is_irreducible(vs::Vector{Contraction}) -> Bool

Check if a set of contractions forms an irreducible graph.
A graph is irreducible if removing any single contraction does not disconnect it.

# Arguments
- `vs`: Vector of contractions representing edges in the graph

# Returns
- `true` if the graph is irreducible (cannot be disconnected by removing one contraction)
- `false` if the graph is reducible (can be disconnected by removing one contraction)

# Algorithm
1. Check if the original graph is connected
2. For each contraction, temporarily remove it from the vector
3. Check if the resulting graph is still connected
4. If any removal disconnects the graph, return false (reducible)
5. If all removals keep the graph connected, return true (irreducible)
"""
function is_irreducible(vs::AbstractVector)
    # If we have fewer than 2 contractions, it's trivially irreducible
    if length(vs) < 2
        return true
    end

    ps = [integer_positions(edge) for edge in vs if is_bulk(edge)] # TODO Keep FixedVector

    # Test removing each contraction
    for i in 1:length(ps)
        # Create a copy without the i-th contraction
        test_contractions = ps[setdiff(begin:end, i)] # TODO SmallCollections.deleteat

        # If removing this contraction disconnects the graph, it's reducible
        if !is_connected(test_contractions)
            return false
        end
    end

    # If we reach here, removing any single contraction keeps the graph connected
    return true
end
