#################################
#       Consistency checks
#################################

is_qq_contraction(v::Contraction) = iszero(sum(Int.(keldysh_index.(v))))
is_qq_contraction(v::Tuple{Field{S},Field{S}}) where {S<:Statistics} =
    is_qq_contraction(Contraction(v))

function is_physical_propagator(a::Contraction)
    positions = position.(a)
    in_out = has_in(positions) ? !has_out(positions) : true
    physical = all(is_physical, a)
    return in_out && physical
end
function is_physical_propagator(a::Tuple{Field{S},Field{S}}) where {S<:Statistics}
    return is_physical_propagator(Contraction(a))
end

function contraction_filter(v::Contraction{S}) where {S<:Statistics}
    out, in = v
    if !contraction_compatible(S, out, in)
        return false
    elseif !is_conserved(v)
        return false
    elseif !is_physical_propagator(v)
        return false
    else
        return !is_qq_contraction(v)
    end
end
function contraction_filter(v::Tuple{Field{S},Field{S}}) where {S<:Statistics}
    return contraction_filter(Contraction(v))
end

"""
    has_two_propagator_zero_loop(vs)

Return whether `vs` contains a zero loop formed by two propagators.
"""
function has_two_propagator_zero_loop(vs::AbstractVector{<:Contraction})
    for i in 1:(length(vs) - 1)
        p1 = positions(vs[i])
        T1 = propagator_type(vs[i]...)
        for j in (i + 1):length(vs)
            p2 = positions(vs[j])
            T2 = propagator_type(vs[j]...)
            if is_retarded(T1) && is_retarded(T2) && is_reversed(p1, p2)
                return true
            elseif is_advanced(T1) && is_advanced(T2) && is_reversed(p1, p2)
                return true
            elseif isequal(p1, p2) && retarded_and_advanced_pair(T1, T2)
                return true
            end
        end
    end
    return false
end

# Position uses Int8: Out is -128, In is 127, and valid bulk positions are
# 1:126. A UInt128 therefore represents every valid position exactly once.
@inline function causal_position_bit(position::Int8)::UInt128
    if position == typemin(Int8)
        return one(UInt128)
    end
    return one(UInt128) << Int(position)
end

@inline function is_causal_propagator(propagator)
    return is_retarded(propagator) || is_advanced(propagator)
end

# Represent a causal propagator as a directed edge source -> destination.
# Gᴿ(x, y) constrains y -> x, while Gᴬ(x, y) constrains x -> y.
@inline function causal_edge(contraction::Contraction, propagator)
    positions = integer_positions(contraction)
    return is_retarded(propagator) ? (positions[2], positions[1]) : positions
end

function visit_causal_position(
    vs::AbstractVector{<:Contraction}, current::Int8, visited::UInt128, active::UInt128
)::Tuple{Bool,UInt128}
    current_bit = causal_position_bit(current)
    visited |= current_bit
    active |= current_bit

    for contraction in vs
        propagator = propagator_type(contraction...)
        is_causal_propagator(propagator) || continue

        source, destination = causal_edge(contraction, propagator)
        source == destination && continue
        source == current || continue

        destination_bit = causal_position_bit(destination)
        !iszero(active & destination_bit) && return true, visited
        if iszero(visited & destination_bit)
            has_loop, visited = visit_causal_position(vs, destination, visited, active)
            has_loop && return true, visited
        end
    end
    return false, visited
end

"""
    has_zero_loop(vs)

Return whether `vs` contains a causal zero loop.
"""
function has_zero_loop(vs::AbstractVector{<:Contraction})
    if length(vs) == 2 && has_two_propagator_zero_loop(vs)
        return true
    end

    visited = zero(UInt128)
    self_positions = zero(UInt128)
    for contraction in vs
        propagator = propagator_type(contraction...)
        is_causal_propagator(propagator) || continue

        source, destination = causal_edge(contraction, propagator)
        if source == destination
            source_bit = causal_position_bit(source)
            if !iszero(self_positions & source_bit)
                return true
            end
            self_positions |= source_bit
            continue
        end

        source_bit = causal_position_bit(source)
        if iszero(visited & source_bit)
            has_loop, visited = visit_causal_position(vs, source, visited, zero(UInt128))
            if has_loop
                return true
            end
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

"""Check whether a propagator satisfies the regularisation convention."""
regular(p::Edge) = regular(Contraction(fields(p)))
regular(p::Tuple{Field{S},Field{S}}) where {S<:Statistics} = regular(Contraction(p))
function regular(qs::Contraction)
    positions = position.(qs)
    if !isequal(positions[1], positions[2])
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

function is_connected(vs::AbstractVector{<:Contraction})
    ps = integer_positions.(vs)
    in_or_out = findfirst(p -> 1 ∈ p || 2 ∈ p, ps) # in case it a vacuum diagram
    edges = isnothing(in_or_out) ? map(p -> p .- 2, ps) : ps
    return is_connected(edges)
end

function is_connected(edges::Union{Vector{Tuple{Int,Int}},Vector{Tuple{Int8,Int8}}})
    all_vertices = vertices(edges)
    if isempty(all_vertices)
        return true
    end

    components = connected_components(all_vertices, edges)
    return length(components) == 1
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
    components = Vector{Set{Int}}()
    remaining = Set(vertices)

    while !isempty(remaining)
        component = Set{Int}()
        queue = [first(remaining)]
        push!(component, first(queue))
        delete!(remaining, first(queue))

        while !isempty(queue)
            current = popfirst!(queue)

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

"""Check whether removing any one bulk contraction disconnects the graph."""
function is_irreducible(vs::AbstractVector)
    if length(vs) < 2
        return true
    end

    ps = [integer_positions(edge) for edge in vs if is_bulk(edge)]

    for i in 1:length(ps)
        test_contractions = ps[setdiff(begin:end, i)]
        if !is_connected(test_contractions)
            return false
        end
    end

    return true
end
