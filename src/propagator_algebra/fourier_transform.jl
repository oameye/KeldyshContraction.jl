"""Numeric coefficient and exact kinematic factor attached to one Fourier graph."""
struct FourierContribution{C<:Number}
    coefficient::C
    kinematic::MomentumPolynomial{ComplexRationals}
end

function Base.isequal(a::FourierContribution{C}, b::FourierContribution{C}) where {C<:Number}
    return isequal(a.coefficient, b.coefficient) && isequal(a.kinematic, b.kinematic)
end
Base.:(==)(a::FourierContribution, b::FourierContribution) = isequal(a, b)
function Base.hash(contribution::FourierContribution, h::UInt)
    return hash(
        FourierContribution,
        hash(contribution.kinematic, hash(contribution.coefficient, h)),
    )
end

"""
Collection of exact momentum-space diagrams.

The dictionary key is a derivative-consumed `FourierDiagram{...,Nothing}` containing
only the canonical coordinate graph and exact routing. The value keeps ordinary numeric
diagram coefficients separate from exact momentum-polynomial kinematics. Different
coordinate derivative placements that lower to the same physical graph therefore meet
under one key without widening the numeric coefficient type.
"""
struct FourierDiagrams{C<:Number,S<:Statistics,E1,E2}
    diagrams::Dict{
        FourierDiagram{S,E1,E2,Nothing},Vector{FourierContribution{C}}
    }
end

function FourierDiagrams{C,S,E1,E2}() where {C<:Number,S<:Statistics,E1,E2}
    K = FourierDiagram{S,E1,E2,Nothing}
    return FourierDiagrams{C,S,E1,E2}(
        Dict{K,Vector{FourierContribution{C}}}()
    )
end

Base.length(collection::FourierDiagrams) = length(collection.diagrams)
Base.isempty(collection::FourierDiagrams) = isempty(collection.diagrams)
Base.iszero(collection::FourierDiagrams) = isempty(collection.diagrams)
Base.iterate(collection::FourierDiagrams) = iterate(collection.diagrams)
Base.iterate(collection::FourierDiagrams, state) = iterate(collection.diagrams, state)
function Base.eltype(::Type{FourierDiagrams{C,S,E1,E2}}) where {C,S,E1,E2}
    return Pair{
        FourierDiagram{S,E1,E2,Nothing},Vector{FourierContribution{C}}
    }
end
function Base.isequal(
    a::FourierDiagrams{C,S,E1,E2}, b::FourierDiagrams{C,S,E1,E2}
) where {C,S,E1,E2}
    return isequal(a.diagrams, b.diagrams)
end
Base.:(==)(a::FourierDiagrams, b::FourierDiagrams) = isequal(a, b)
Base.hash(collection::FourierDiagrams, h::UInt) = hash(collection.diagrams, h)

function _push_fourier!(
    collection::FourierDiagrams{C,S,E1,E2},
    graph::FourierDiagram{S,E1,E2,Nothing},
    coefficient::Number,
    kinematic::MomentumPolynomial{ComplexRationals},
) where {C<:Number,S<:Statistics,E1,E2}
    value = _simplify(convert(C, coefficient))
    iszero(value) && return collection

    contributions = get!(collection.diagrams, graph) do
        FourierContribution{C}[]
    end
    for i in eachindex(contributions)
        existing = contributions[i]
        isequal(existing.kinematic, kinematic) || continue
        combined = _simplify(existing.coefficient + value)
        if iszero(combined)
            deleteat!(contributions, i)
            isempty(contributions) && delete!(collection.diagrams, graph)
        else
            contributions[i] = FourierContribution{C}(combined, kinematic)
        end
        return collection
    end
    push!(contributions, FourierContribution{C}(value, kinematic))
    return collection
end

@inline function _without_derivatives(field::Field{S}) where {S<:Statistics}
    return reconstruct(field; derivative=DerivativeMultiIndex())
end

function _without_derivatives(edge::Edge{S}) where {S<:Statistics}
    return Edge{S}(
        _without_derivatives(edge.out),
        _without_derivatives(edge.in),
        edge.edgetype,
        Momenta(),
    )
end

@inline function _field_fourier_key(field::Field)
    return (
        name(field),
        slots(field_indices(field)),
        Int(orientation(field)),
        Int(keldysh_index(field)),
        Int(regularisation(field)),
    )
end

@inline function _edge_fourier_key(edge::Edge)
    return (
        sort_by_position_and_type(edge),
        _field_fourier_key(edge.out),
        _field_fourier_key(edge.in),
        Int(propagator_type(edge)),
    )
end

@inline function _derivative_tiebreak(edge::Edge)
    return (
        derivative_multiindex(edge.out).orders,
        derivative_multiindex(edge.in).orders,
    )
end

function _fourier_source_isless(
    stripped_a::Edge,
    decorated_a::Edge,
    stripped_b::Edge,
    decorated_b::Edge,
)
    key_a = _edge_fourier_key(stripped_a)
    key_b = _edge_fourier_key(stripped_b)
    key_a == key_b || return isless(key_a, key_b)
    derivative_a = _derivative_tiebreak(decorated_a)
    derivative_b = _derivative_tiebreak(decorated_b)
    derivative_a == derivative_b && return false
    return isless(derivative_a, derivative_b)
end

function _canonical_fourier_source(
    diagram::Diagram{S,E1,E2},
) where {S<:Statistics,E1,E2}
    decorated = Edge{S}[
        Edge{S}(edge.out, edge.in, edge.edgetype, Momenta()) for
        edge in contractions(diagram)
    ]
    stripped = Edge{S}[_without_derivatives(edge) for edge in decorated]

    graph_positions = canonicalization_positions(stripped)
    physical_permutation, _, _, _ =
        canonicalization_permutations(stripped, graph_positions)
    mapping = make_permutation_dict(physical_permutation, graph_positions, stripped)

    decorated = Edge{S}[relabel_bulk_positions(edge, mapping) for edge in decorated]
    stripped = Edge{S}[relabel_bulk_positions(edge, mapping) for edge in stripped]

    permutation = sortperm(
        collect(eachindex(stripped));
        lt=(i, j) -> _fourier_source_isless(
            stripped[i], decorated[i], stripped[j], decorated[j]
        ),
    )
    decorated = decorated[permutation]
    stripped = stripped[permutation]

    stripped_fixed = FixedVector{E1,Edge{S}}(stripped)
    stripped_diagram = Diagram(stripped_fixed, Val(E2))
    decorated_fixed = FixedVector{E1,Edge{S}}(decorated)
    decorated_diagram = Diagram{S,E1,E2}(decorated_fixed, topology(stripped_diagram))

    routed = FourierDiagram(stripped_diagram)
    decorated_routed = FourierDiagram{S,E1,E2,Nothing}(
        decorated_diagram,
        routed.basis,
        routed.edge_momenta,
        routed.external_count,
        routed.loop_count,
        nothing,
    )
    lowered = lower_fourier_derivatives(decorated_routed)
    return routed, kinematic_factor(lowered), lowered
end

"""Transform one coordinate-space diagram to its exact routed/lowered Fourier form."""
function fourier_transform(diagram::Diagram{S,E1,E2}) where {S<:Statistics,E1,E2}
    _, _, lowered = _canonical_fourier_source(diagram)
    return lowered
end

"""Transform and derivative-consume a collection of coordinate-space diagrams."""
function fourier_transform(
    diagrams::Diagrams{C,S,E1,E2},
) where {C<:Number,S<:Statistics,E1,E2}
    out = FourierDiagrams{C,S,E1,E2}()
    for (diagram, coefficient) in diagrams
        graph, kinematic, _ = _canonical_fourier_source(diagram)
        _push_fourier!(out, graph, coefficient, kinematic)
    end
    return out
end

"""Exact Fourier-space dressed propagator."""
struct FourierDressedPropagator{C<:Number,S<:Statistics,O,E1,E2}
    keldysh::FourierDiagrams{C,S,E1,E2}
    retarded::FourierDiagrams{C,S,E1,E2}
    advanced::FourierDiagrams{C,S,E1,E2}
    parameter::ParameterMonomial
end

order(::FourierDressedPropagator{C,S,O}) where {C,S,O} = O
statistics(::FourierDressedPropagator{C,S}) where {C,S} = S
parameters(G::FourierDressedPropagator) = G.parameter
function Base.isequal(
    a::FourierDressedPropagator{C,S,O,E1,E2},
    b::FourierDressedPropagator{C,S,O,E1,E2},
) where {C,S,O,E1,E2}
    return isequal(a.keldysh, b.keldysh) &&
           isequal(a.retarded, b.retarded) &&
           isequal(a.advanced, b.advanced) &&
           isequal(a.parameter, b.parameter)
end
Base.:(==)(a::FourierDressedPropagator, b::FourierDressedPropagator) = isequal(a, b)
Base.hash(G::FourierDressedPropagator, h::UInt) =
    hash((G.keldysh, G.retarded, G.advanced, G.parameter), h)

function fourier_transform(
    G::DressedPropagator{C,S,O,E1,E2},
) where {C<:Number,S<:Statistics,O,E1,E2}
    return FourierDressedPropagator{C,S,O,E1,E2}(
        fourier_transform(G.keldysh),
        fourier_transform(G.retarded),
        fourier_transform(G.advanced),
        G.parameter,
    )
end

"""Exact Fourier-space self-energy."""
struct FourierSelfEnergy{C<:Number,S<:Statistics,O,E1,E2}
    keldysh::FourierDiagrams{C,S,E1,E2}
    retarded::FourierDiagrams{C,S,E1,E2}
    advanced::FourierDiagrams{C,S,E1,E2}
    parameter::ParameterMonomial
end

order(::FourierSelfEnergy{C,S,O}) where {C,S,O} = O
statistics(::FourierSelfEnergy{C,S}) where {C,S} = S
parameters(Σ::FourierSelfEnergy) = Σ.parameter
function Base.isequal(
    a::FourierSelfEnergy{C,S,O,E1,E2},
    b::FourierSelfEnergy{C,S,O,E1,E2},
) where {C,S,O,E1,E2}
    return isequal(a.keldysh, b.keldysh) &&
           isequal(a.retarded, b.retarded) &&
           isequal(a.advanced, b.advanced) &&
           isequal(a.parameter, b.parameter)
end
Base.:(==)(a::FourierSelfEnergy, b::FourierSelfEnergy) = isequal(a, b)
Base.hash(Σ::FourierSelfEnergy, h::UInt) =
    hash((Σ.keldysh, Σ.retarded, Σ.advanced, Σ.parameter), h)

function _amputate_fourier_graph(
    graph::FourierDiagram{S,E,E2,Nothing}, ::Val{SE}, ::Val{ST}
) where {S<:Statistics,E,E2,SE,ST}
    edges = Edge{S}[]
    momenta = LinearMomentum[]
    for (edge, momentum) in zip(contractions(graph.coordinate), graph.edge_momenta)
        is_bulk(edge) || continue
        push!(edges, edge)
        push!(momenta, momentum)
    end
    length(edges) == SE || error("Fourier self-energy amputation found the wrong bulk-edge count")

    fixed_edges = FixedVector{SE,Edge{S}}(edges)
    coordinate = Diagram(fixed_edges, Val(ST))
    fixed_momenta = FixedVector{SE,LinearMomentum}(momenta)
    return FourierDiagram{S,SE,ST,Nothing}(
        coordinate,
        graph.basis,
        fixed_momenta,
        graph.external_count,
        graph.loop_count,
        nothing,
    )
end

function construct_fourier_self_energy!(
    self_energy::SmallCollections.SmallDict,
    diagrams::FourierDiagrams{C,S,E,E2},
    ::Val{SE},
    ::Val{ST},
) where {C<:Number,S<:Statistics,E,E2,SE,ST}
    for (graph, contributions) in diagrams
        coordinate_edges = contractions(graph.coordinate)
        is_irreducible(coordinate_edges) || continue

        categories = position_category.(coordinate_edges)
        types = propagator_type.(coordinate_edges)
        external_types = SmallCollections.SmallDict{E,Symbol,PropagatorType.T}(
            category => propagator for (category, propagator) in zip(categories, types)
        )
        if is_keldysh(external_types[:out]) && is_keldysh(external_types[:in])
            continue
        end

        component = self_energy_type(S, external_types)
        amputated = _amputate_fourier_graph(graph, Val(SE), Val(ST))
        for contribution in contributions
            _push_fourier!(
                self_energy[component],
                amputated,
                contribution.coefficient,
                contribution.kinematic,
            )
        end
    end
    return self_energy
end

function _fourier_self_energy(
    G::FourierDressedPropagator{C,S,O,E1,E2},
) where {C<:Number,S<:Statistics,O,E1,E2}
    SE = E1 - 2
    ST = max_edges(O)
    D = FourierDiagrams{C,S,SE,ST}
    self_energy = SmallCollections.SmallDict{3,PropagatorType.T,D}((
        PropagatorType.Advanced => D(),
        PropagatorType.Retarded => D(),
        PropagatorType.Keldysh => D(),
    ))
    construct_fourier_self_energy!(self_energy, G.keldysh, Val(SE), Val(ST))
    return FourierSelfEnergy{C,S,O,SE,ST}(
        self_energy[PropagatorType.Keldysh],
        self_energy[PropagatorType.Retarded],
        self_energy[PropagatorType.Advanced],
        G.parameter,
    )
end

SelfEnergy(G::FourierDressedPropagator) = _fourier_self_energy(G)

@inline function _diagram_has_derivatives(diagram::Diagram)
    for edge in contractions(diagram), field in fields(edge)
        isempty(derivative_multiindex(field)) || return true
    end
    return false
end

function _collection_has_derivatives(diagrams::Diagrams)
    return any(_diagram_has_derivatives(diagram) for diagram in keys(diagrams.diagrams))
end

"""
Transform a coordinate-space self-energy when its amputated representation is derivative-free.

For derivative-decorated interactions, coordinate-space amputation can already have removed
external endpoint derivative metadata. In that case transform the `DressedPropagator` first
and construct `SelfEnergy(fourier_transform(G))` so the external momentum factors are retained.
"""
function fourier_transform(
    Σ::SelfEnergy{C,S,O,E1,E2},
) where {C<:Number,S<:Statistics,O,E1,E2}
    if _collection_has_derivatives(Σ.keldysh) ||
       _collection_has_derivatives(Σ.retarded) ||
       _collection_has_derivatives(Σ.advanced)
        throw(
            ArgumentError(
                "cannot safely Fourier-transform a derivative-decorated coordinate SelfEnergy after amputation; transform the DressedPropagator first and construct SelfEnergy(fourier_transform(G))",
            ),
        )
    end
    return FourierSelfEnergy{C,S,O,E1,E2}(
        fourier_transform(Σ.keldysh),
        fourier_transform(Σ.retarded),
        fourier_transform(Σ.advanced),
        Σ.parameter,
    )
end

function _fourier_zero(::Type{FourierDiagrams{C,S,E1,E2}}) where {C,S,E1,E2}
    return FourierDiagrams{C,S,E1,E2}()
end

function matrix(
    G::FourierDressedPropagator{C,Boson,O,E1,E2},
) where {C<:Number,O,E1,E2}
    D = FourierDiagrams{C,Boson,E1,E2}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = G.keldysh
    result[1, 2] = G.retarded
    result[2, 1] = G.advanced
    result[2, 2] = _fourier_zero(D)
    return result
end

function matrix(
    G::FourierDressedPropagator{C,Fermion,O,E1,E2},
) where {C<:Number,O,E1,E2}
    D = FourierDiagrams{C,Fermion,E1,E2}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = G.retarded
    result[1, 2] = G.keldysh
    result[2, 1] = _fourier_zero(D)
    result[2, 2] = G.advanced
    return result
end

function matrix(
    Σ::FourierSelfEnergy{C,Boson,O,E1,E2},
) where {C<:Number,O,E1,E2}
    D = FourierDiagrams{C,Boson,E1,E2}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = _fourier_zero(D)
    result[1, 2] = Σ.advanced
    result[2, 1] = Σ.retarded
    result[2, 2] = Σ.keldysh
    return result
end

function matrix(
    Σ::FourierSelfEnergy{C,Fermion,O,E1,E2},
) where {C<:Number,O,E1,E2}
    D = FourierDiagrams{C,Fermion,E1,E2}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = Σ.retarded
    result[1, 2] = Σ.keldysh
    result[2, 1] = _fourier_zero(D)
    result[2, 2] = Σ.advanced
    return result
end
