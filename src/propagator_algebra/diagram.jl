struct Diagram{E1,E2}
    contractions::FixedVector{E1,Edge}
    topology::FixedVector{E2,Int}
end

function Diagram(
    contractions::Vector{T}, ::Val{E}, ::Val{E2}
) where {T<:Union{Contraction,Edge},E,E2}
    @assert length(contractions) == E "The supplied Val{edges} must match the contraction count"
    @assert E > 0 "Contraction vector must not be empty"
    sort!(contractions; by=sort_by_position_and_type)

    edges = SmallCollections.FixedVector{E,Edge}(
        T <: Contraction ? Edge(c) : c for c in contractions
    )
    return Diagram(edges, Val(E2))
end

function Diagram(edges::FixedVector{E,Edge}, ::Val{E2}) where {E,E2}
    topology = bulk_multiplicity(edges, Val(E2))
    @assert length(topology) == E2 "The supplied Val{topology} must match the topology size"
    topology = SmallCollections.FixedVector{E2,Int}(topology)
    return Diagram{E,E2}(edges, topology)
end

function Diagram(d::Diagram{E1,E2}, momenta::FixedVector{E1,Momenta}) where {E1,E2}
    contractions = FixedVector{E1,Edge}(
        Edge(c, m) for (c, m) in zip(d.contractions, momenta)
    )
    return Diagram{E1,E2}(contractions, d.topology)
end

Base.isequal(d1::Diagram, d2::Diagram) = isequal(contractions(d1), contractions(d2))
Base.hash(d::Diagram, h::UInt) = hash(contractions(d), h)
contractions(d::Diagram) = d.contractions
topology(d::Diagram) = d.topology
function momenta(d::Diagram)
    return map(momenta, d.contractions)
end
topology_length(x::Int) = max(0, x - 4) # TODO review this

Base.isempty(d::Diagram) = isempty(d.contractions)

function set_reg_to_zero(d::Diagram{E1,E2}) where {E1,E2}
    new_contractions = map(set_reg_to_zero, d.contractions)
    return Diagram{E1,E2}(new_contractions, d.topology)
end

################
#   Diagrams
###############

struct Diagrams{E1,E2}
    diagrams::Dict{Diagram{E1,E2},ComplexRationals}
end # TODO try SwissDict or RobinDict from DataStructures.jl.
function Diagrams{E1,E2}() where {E1,E2}
    dict = Dict{Diagram{E1,E2},ComplexRationals}()
    return Diagrams(dict)
end
function Diagrams(
    diagrams::Vector{Diagram{E1,E2}}, prefactor::ComplexRationals
) where {E1,E2}
    dict = Dict{Diagram{E1,E2},ComplexRationals}(d => prefactor for d in diagrams)
    return Diagrams{E1,E2}(dict)
end
function Diagrams(
    contractions::Vector{Vector{Contraction}},
    prefactor::ComplexRationals,
    ::Val{E},
    ::Val{E2},
) where {E,E2}
    @assert length(contractions) > 0 "Contraction vector must not be empty"
    c = first(contractions)
    @assert length(c) == E "The supplied Val{edges} must match the contraction count"

    imag_factor = im^E # Contraction becomes propagator
    dict = Dict{Diagram{E,E2},ComplexRationals}(
        Diagram(c, Val(E), Val(E2)) => _simplify(imag_factor * prefactor) for
        c in contractions
    )
    return Diagrams{E,E2}(dict)
end
Base.isequal(d1::Diagrams, d2::Diagrams) = isequal(d1.diagrams, d2.diagrams)
Base.hash(d::Diagrams, h::UInt) = hash(d.diagrams, h)
Base.iszero(d::Diagrams) = isempty(d.diagrams)
SmallCollections.default(::Type{Diagrams}) = Diagrams{0,0}()
SmallCollections.default(::Type{Diagrams{E1,E2}}) where {E1,E2} = Diagrams{E1,E2}()

number_of_propagators(a::QMul) = length(a) ÷ 2
number_of_propagators(a::QAdd) = length(first(a.arguments)) ÷ 2
number_of_propagators(L::InteractionLagrangian) = length(first(L.lagrangian.arguments)) ÷ 2

# Add a single diagram, summing prefactors if it already exists
function Base.push!(collection::Diagrams, diagram::Diagram, prefactor::Number)
    if haskey(collection.diagrams, diagram)
        collection.diagrams[diagram] += prefactor
        if iszero(collection.diagrams[diagram])
            delete!(collection.diagrams, diagram)
        end
    else
        collection.diagrams[diagram] = prefactor
    end
    return collection
end

function filter_nonzero!(collection::Diagrams)
    filter!((kv) -> !iszero(kv[2]), collection.diagrams)
    return collection
end

# Convert to vector of diagrams (ignoring prefactors)
function Base.collect(collection::Diagrams)
    return collect(keys(collection.diagrams))
end

function Base.:*(prefactor::Number, collection::Diagrams)
    collection = deepcopy(collection) # needed so that -1 * Σ.keldysh does not change Σ
    foreach(collection) do (diagram, _)
        collection.diagrams[diagram] *= prefactor
        return collection.diagrams[diagram] = _simplify(collection.diagrams[diagram])
    end
    return collection
end
Base.:*(diagrams::Diagrams, prefactor::Number) = prefactor * diagrams

# Make the collection iterable (iterate over pairs)
Base.iterate(collection::Diagrams) = iterate(collection.diagrams)
Base.iterate(collection::Diagrams, state) = iterate(collection.diagrams, state)
Base.length(collection::Diagrams) = length(collection.diagrams)
Base.eltype(::Type{Diagrams}) = Pair{Diagram,Number}
Base.eltype(::Type{Diagrams{E1,E2}}) where {E1,E2} = Pair{Diagram{E1,E2},ComplexRationals}

function Base.adjoint(d::Diagrams)
    dict = Dict(adjoint_diagram(pair) for pair in d)
    return Diagrams(dict)
end
"""
    has_external_legs(d::Diagram)

Whether the diagram carries [`In`](@ref) and [`Out`](@ref) coordinates. A
[`DressedPropagator`](@ref) diagram does; an amputated [`SelfEnergy`](@ref) diagram does
not, so there is no way to tell which of its bulk coordinates the external legs attach to.
"""
has_external_legs(d::Diagram) = any(e -> is_in(fields(e)), contractions(d))

"""
    amputated_leg_indices(d::Diagram)

The bulk coordinates an amputated [`SelfEnergy`](@ref) diagram lost its external legs to,
as `(out_leg, in_leg)`, or `(0, 0)` when there are none to find. Every interaction vertex
of an [`InteractionLagrangian`](@ref) has the same valence, so the vertex one [`Create`](@ref)
short is where the outgoing leg was cut and the vertex one [`Destroy`](@ref) short is where
the incoming one was. Both legs may hang off the same vertex.
"""
function amputated_leg_indices(d::Diagram)
    _contractions = contractions(d)
    verts = sort!(unique(Int(index(f)) for e in _contractions for f in fields(e)))
    valence, remainder = divrem(length(_contractions) + 1, length(verts))
    iszero(remainder) || return (0, 0)
    outgoing = zeros(Int, length(verts))
    incoming = zeros(Int, length(verts))
    for e in _contractions
        _out, _in = fields(e)
        outgoing[searchsortedfirst(verts, Int(index(_out)))] += 1
        incoming[searchsortedfirst(verts, Int(index(_in)))] += 1
    end
    out_leg = count(<(valence), incoming) == 1 ? verts[findfirst(<(valence), incoming)] : 0
    in_leg = count(<(valence), outgoing) == 1 ? verts[findfirst(<(valence), outgoing)] : 0
    return (out_leg, in_leg)
end

"""
    with_external_legs(d::Diagram)

The contractions of `d` together with the external legs it is missing, and how many were
added. Reversing a diagram exchanges its incoming and outgoing coordinate, so an amputated
diagram has to get its legs back before it can be reversed at all.
"""
function with_external_legs(d::Diagram)
    _contractions = Contraction[fields(e) for e in contractions(d)]
    has_external_legs(d) && return (_contractions, 0)
    out_leg, in_leg = amputated_leg_indices(d)
    (iszero(out_leg) || iszero(in_leg)) && return (_contractions, -1)
    _out, _in = first(_contractions)
    push!(_contractions, (_out(Out()), _in(Bulk(out_leg))))
    push!(_contractions, (_out(Bulk(in_leg)), _in(In())))
    return (_contractions, 2)
end

function adjoint_diagram(
    pair::Pair{Diagram{E1,E2},ComplexRationals}
)::Pair{Diagram{E1,E2},ComplexRationals} where {E1,E2}
    # G^R = G^A
    # G^K = -G^K
    d, prefactor = pair
    _contractions = contractions(d)
    minus_signs = count(is_keldysh, _contractions)
    prefactor′ = prefactor * (-1)^minus_signs
    legged, added = with_external_legs(d)
    adjoint_edges = if added < 0
        # Nothing to reverse around: fall back to swapping the contours in place.
        Edge[adjoint(e) for e in _contractions]
    else
        # Reversing every edge exchanges the two external coordinates, so relabel the bulk
        # coordinates back into canonical form, then drop the legs again.
        reversed = Contraction[reverse_contraction(c) for c in legged]
        canonical = canonicalize(reversed)
        kept = Contraction[canonical[i] for i in 1:(length(canonical) - added)]
        Edge[Edge(c) for c in kept]
    end
    sorted_adjoint_edges = sort(collect(adjoint_edges); by=sort_by_position_and_type)
    edges = SmallCollections.FixedVector{E1,Edge}(e for e in sorted_adjoint_edges)
    return Diagram(edges, Val(E2)) => _simplify(adjoint(prefactor′))
end

function set_reg_to_zero(d::Diagrams{E1,E2}) where {E1,E2}
    diagrams = Diagrams{E1,E2}()
    for (diagram, value) in d
        push!(diagrams, set_reg_to_zero(diagram), value)
    end
    return diagrams
end

"""
Compute multiplicity of the edges between two different vertices in the bulk.
The resulting vector is a signature for the topology of the diagram.
"""
function bulk_multiplicity(edges::AbstractArray{Tuple{Int8,Int8}})
    edges = filter(is_not_equal_time_bulk_edge, edges)

    vert = vertices(edges)
    m = max_edges(length(vert))

    mult = zeros(Int, m)
    for edge in edges
        idx = edge_to_index(edge[1], edge[2], length(vert))
        mult[idx] += 1
    end
    return mult
end
function bulk_multiplicity(vs::AbstractArray{Edge})
    return bulk_multiplicity(map(integer_positions, vs))
end

function bulk_multiplicity(edges::AbstractArray{Tuple{Int8,Int8}}, ::Val{E2}) where {E2}
    edges = filter(is_not_equal_time_bulk_edge, edges)

    vert = vertices(edges)
    m = max_edges(length(vert))
    @assert m == E2 "The supplied Val{topology} must match the topology size"

    if iszero(E2)
        mult = SmallCollections.MutableFixedVector{0,Int}(undef)
    else
        mult = SmallCollections.MutableFixedVector{E2,Int}(0 for _ in 1:E2)
    end
    for edge in edges
        idx = edge_to_index(edge[1], edge[2], length(vert))
        mult[idx] += 1
    end
    return SmallCollections.FixedVector(mult)
end

function bulk_multiplicity(vs::AbstractArray{Edge}, ::Val{E2}) where {E2}
    return bulk_multiplicity(map(integer_positions, vs), Val(E2))
end

max_edges(n::Integer)::Integer = n * (n - 1) ÷ 2

function is_not_equal_time_bulk_edge(edge)
    return !(typemin(Int8) ∈ edge) && !(typemax(Int8) ∈ edge) && !isequal(edge[1], edge[2])
end

function topologies(ds::Diagrams{E1,E2}) where {E1,E2}
    terms = collect(keys(ds.diagrams))
    diagram_topologies = getfield.(terms, :topology)
    _topologies = unique(diagram_topologies)
    topologies = Dict{FixedVector{E2,Int},Vector{Diagram{E1,E2}}}()
    for topology in _topologies
        idxs = findall(
            i ->
                length(i) == length(topology) &&
                all(j -> i[j] == topology[j], eachindex(i)),
            diagram_topologies,
        )
        topologies[topology] = terms[idxs]
    end
    return topologies
end

function _simplify_prefactors!(g::Diagrams)
    for k in keys(g.diagrams)
        g.diagrams[k] = _simplify(g.diagrams[k])
    end
    return nothing
end
