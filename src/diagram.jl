struct Diagram{E1,E2}
    contractions::SmallCollections.FixedVector{E1,Edge}
    topology::SmallCollections.FixedVector{E2,Int64}
end

function Diagram(contractions::Vector{T}) where {T<:Union{Contraction,Edge}}
    @assert length(contractions) > 0 "Contraction vector must not be empty"
    sort!(contractions; by=sort_by_position_and_type)

    edges = SmallCollections.FixedVector{length(contractions),Edge}(
        T <: Contraction ? Edge(c) : c for c in contractions
    )
    # TODO this is type unstable move to FixedSizeArrays.jl when released
    # https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/115
    return Diagram(edges)
end

function Diagram(edges::FixedVector{E,Edge}) where {E}
    topology = bulk_multiplicity(edges)
    Diagram{E,length(topology)}(edges, topology)
end

Base.isequal(d1::Diagram, d2::Diagram) = isequal(contractions(d1), contractions(d2))
Base.hash(d::Diagram, h::UInt) = hash(contractions(d), h)
contractions(d::Diagram) = d.contractions
topology(d::Diagram) = d.topology

topology_length(x::Int) = max(0, x - 4)

################
#   Diagrams
###############

struct Diagrams{E1,E2}
    diagrams::Dict{Diagram{E1,E2},ComplexF64}
end # TODO try SwissDict or RobinDict from DataStructures.jl.
function Diagrams{E1,E2}() where {E1,E2}
    dict = Dict{Diagram{E1,E2},ComplexF64}()
    return Diagrams(dict)
end
function Diagrams(diagrams::Vector{Diagram{E1,E2}}, prefactor::ComplexF64) where {E1,E2}
    dict = Dict{Diagram{E1,E2},ComplexF64}(d => prefactor for d in diagrams)
    return Diagrams{E1,E2}(dict)
end
function Diagrams(contractions::Vector{Vector{Contraction}}, prefactor::ComplexF64)
    @assert length(contractions) > 0 "Contraction vector must not be empty"
    c = first(contractions)
    E = length(c)

    imag_factor = im^E # Contraction becomes propagator
    dict = Dict{Diagram{E,topology_length(E)},ComplexF64}(
        Diagram(c) => _simplify(imag_factor * prefactor) for c in contractions
    )
    return Diagrams{E,topology_length(E)}(dict)
end
Base.isequal(d1::Diagrams, d2::Diagrams) = isequal(d1.diagrams, d2.diagrams)
Base.hash(d::Diagrams, h::UInt) = hash(d.diagrams, h)
Base.iszero(d::Diagrams) = isempty(d.diagrams)
SmallCollections.default(::Type{Diagrams}) = Diagrams{0,0}()

number_of_propagators(a::QMul) = length(a) ÷ 2
number_of_propagators(L::InteractionLagrangian) = length(first(L.lagrangian.arguments)) ÷ 2

# Add a single diagram, summing prefactors if it already exists
function Base.push!(collection::Diagrams, diagram::Diagram, prefactor::Number)
    if haskey(collection.diagrams, diagram)
        collection.diagrams[diagram] += prefactor
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
function multiply!(collection::Diagrams, prefactor::Number)
    foreach(collection) do (diagram, coeff)
        collection.diagrams[diagram] *= prefactor
        collection.diagrams[diagram] = _simplify(collection.diagrams[diagram])
    end
    return collection
end # TODO Instead of multiply! we could use Base.:*
function Base.:*(prefactor::Number, collection::Diagrams)
    diagrams = deepcopy(collection)
    multiply!(diagrams, prefactor)
    return diagrams
end

# Make the collection iterable (iterate over pairs)
Base.iterate(collection::Diagrams) = iterate(collection.diagrams)
Base.iterate(collection::Diagrams, state) = iterate(collection.diagrams, state)
Base.length(collection::Diagrams) = length(collection.diagrams)
Base.eltype(::Type{Diagrams}) = Pair{Diagram,Number}

function Base.adjoint(d::Diagrams)
    dict = Dict(adjoint_diagram(pair) for pair in d)
    return Diagrams(dict)
end
function adjoint_diagram(
    pair::Pair{Diagram{E1,E2},ComplexF64}
)::Pair{Diagram{E1,E2},ComplexF64} where {E1,E2}
    # G^R = G^A
    # G^K = -G^K
    d, prefactor = pair
    _contractions = contractions(d)
    minus_signs = count(is_keldysh, _contractions)
    prefactor′ = prefactor * (-1)^minus_signs
    adjoint_edges = Edge[adjoint(e) for e in _contractions]
    sorted_adjoint_edges = sort(collect(adjoint_edges); by=sort_by_position_and_type)
    edges = SmallCollections.FixedVector{E1,Edge}(e for e in sorted_adjoint_edges)
    return Diagram(edges) => _simplify(adjoint(prefactor′))
end

"""
Compute multiplicity of the edges between two different vertices in the bulk.
The resulting vector is a signature for the topology of the diagram.
"""
function bulk_multiplicity(edges::FixedVector{N,Tuple{Int8,Int8}}) where {N}
    edges = filter(is_not_equal_time_bulk_edge, edges)

    vert = vertices(edges)
    m = max_edges(length(vert))

    if iszero(m)
        mult = SmallCollections.MutableFixedVector{0,Int}(undef)
    else
        mult = SmallCollections.MutableFixedVector{m,Int}(0 for i in 1:m)
    end # https://github.com/matthias314/SmallCollections.jl/issues/12
    for edge in edges
        idx = edge_to_index(edge[1], edge[2], length(vert))
        mult[idx] += 1
    end
    return SmallCollections.FixedVector(mult)
end
function bulk_multiplicity(vs::AbstractArray{Edge})
    return bulk_multiplicity(map(integer_positions, vs))
end

max_edges(n::Integer)::Integer = n * (n - 1) ÷ 2

function is_not_equal_time_bulk_edge(edge)
    return !(typemin(Int8) ∈ edge) && !(typemax(Int8) ∈ edge) && !isequal(edge[1], edge[2])
end

function topologies(ds::Diagrams)
    terms = collect(keys(ds.diagrams))
    _bulk_multiplicity = map(terms) do diagram
        vs = map(diagram.contractions) do c
            ff = fields(c)
            integer_positions((ff[1], ff[2]))
        end
        bulk_multiplicity(vs)
    end
    _topologies = unique(_bulk_multiplicity)
    pairs = map(_topologies) do t
        idxs = findall(i -> i == t, _bulk_multiplicity)
        t => terms[idxs]
    end
    Dict(pairs)
end

function _simplify_prefactors!(g::Diagrams)
    for k in keys(g.diagrams)
        g.diagrams[k] = _simplify(g.diagrams[k])
    end
    return nothing
end
