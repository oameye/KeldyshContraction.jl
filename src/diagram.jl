struct Diagram{E,T}
    contractions::T
    topology::Vector{Int64} # TODO: make immutable
end

function Diagram(contractions::Vector{<:Contraction})
    @assert length(contractions) > 0 "Contraction vector must not be empty"
    E = length(contractions)
    sort!(contractions; by=sort_by_position_and_type)
    edges = StaticArrays.sacollect(
        SVector{length(contractions),Edge}, Edge(c) for c in contractions
    ) # TODO this is type unstable move to FixedSizeArrays.jl when released
    # https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/115
    return Diagram(edges)
end
function Diagram(contractions::Vector{<:Edge})
    @assert length(contractions) > 0 "Contraction vector must not be empty"
    E = length(contractions)
    sort!(contractions; by=sort_by_position_and_type)
    # TODO: sort to be sure?
    edges = StaticArrays.sacollect(
        SVector{length(contractions),Edge}, c for c in contractions
    )
    return Diagram(edges)
end # TODO: make above two functions less redundant
function Diagram(edges::SVector{E,Edge}) where {E}
    Diagram{E,SVector{E,Edge}}(edges, bulk_multiplicity(edges))
end
Base.isequal(d1::Diagram, d2::Diagram) = isequal(contractions(d1), contractions(d2))
Base.hash(d::Diagram, h::UInt) = hash(contractions(d), h)
contractions(d::Diagram) = d.contractions
topology(d::Diagram) = d.topology

################
#   Diagrams
###############

struct Diagrams{E,T}
    diagrams::Dict{Diagram{E,T},ComplexF64}
end # TODO try SwissDict or RobinDict from DataStructures.jl.
function Diagrams{E}() where {E}
    dict = Dict{Diagram{E,SVector{E,Edge}},ComplexF64}()
    return Diagrams(dict)
end
function Diagrams(diagrams::Vector{Diagram{E,T}}, prefactor::ComplexF64) where {E,T}
    dict = Dict{Diagram{E,T},ComplexF64}(d => prefactor for d in diagrams)
    return Diagrams{E,T}(dict)
end
function Diagrams(contractions::Vector{Vector{Contraction}}, prefactor::ComplexF64)
    @assert length(contractions) > 0 "Contraction vector must not be empty"
    c = first(contractions)
    E = length(c)
    T = SVector{E,Edge}
    imag_factor = im^E # Contraction becomes propagator
    dict = Dict{Diagram{E,T},ComplexF64}(
        Diagram(c) => _simplify(imag_factor * prefactor) for c in contractions
    )
    return Diagrams{E,T}(dict)
end
Base.isequal(d1::Diagrams, d2::Diagrams) = isequal(d1.diagrams, d2.diagrams)
Base.hash(d::Diagrams, h::UInt) = hash(d.diagrams, h)
Base.iszero(d::Diagrams) = isempty(d.diagrams)

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
end # TODO Instead of push! we could use setindex! to be more consistent with the rest of Julia

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

function Base.adjoint(d::Diagrams{E}) where {E}
    dict = Dict(adjoint_diagram(pair) for pair in d)
    return Diagrams(dict)
end
function adjoint_diagram(
    pair::Pair{Diagram{E,T},ComplexF64}
)::Pair{Diagram{E,T},ComplexF64} where {E,T}
    # G^R = G^A
    # G^K = -G^K
    d, prefactor = pair
    _contractions = contractions(d)
    minus_signs = count(is_keldysh, _contractions)
    prefactor′ = prefactor * (-1)^minus_signs
    adjoint_edges = Edge[adjoint(e) for e in _contractions]
    sorted_adjoint_edges = sort(collect(adjoint_edges); by=sort_by_position_and_type)
    edges = StaticArrays.sacollect(T, e for e in sorted_adjoint_edges)
    return Diagram(edges) => _simplify(adjoint(prefactor′))
end

function bulk_multiplicity(edges::AbstractArray{Tuple{Int,Int}})
    ff = edge -> !(1 ∈ edge) && !(2 ∈ edge) && !isequal(edge[1], edge[2])
    edges = filter(ff, edges)
    map!(edge -> edge .- 2, edges, edges)

    vert = vertices(edges)
    m = max_edges(length(vert))
    mult = SmallCollections.MutableSmallVector{m,Int}(0 for i in 1:m)
    for edge in edges # TODO: probably should use pairing function
        idx = edge_to_index(edge[1], edge[2], length(vert))
        mult[idx] += 1
    end
    return mult
end # TODO: do we need to sort?
function bulk_multiplicity(vs::AbstractArray{Edge})
    return bulk_multiplicity(map(integer_positions, vs))
end

max_edges(n::Int)::Int = n * (n - 1) ÷ 2

"""
maps edge (i,j) to a unique integer in range 1:max_edges(n).
It assumes that i ≠ j
"""
function edge_to_index(i::Int, j::Int, n::Int)
    # Ensure i < j
    i, j = minmax(i, j)

    # Calculate unique index
    # This maps edge (i,j) to a unique integer in range 1:max_edges(n)
    return (i - 1) * (n - i ÷ 2) + (j - i)
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
