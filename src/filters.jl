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
function has_zero_loop(vs::Vector{Contraction})
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

function is_reversed(p1, p2)
    return isequal(p1[1], p2[2]) && isequal(p1[2], p2[1])
end
function retarded_and_advanced_pair(T1, T2)
    return (is_retarded(T1) && is_advanced(T2)) || (is_advanced(T1) && is_retarded(T2))
end
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
    elseif subtraction(_reg) < 0 && T == Retarded
        return false
    elseif subtraction(_reg) > 0 && T == Advanced
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
