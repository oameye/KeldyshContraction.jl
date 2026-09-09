"""Numeric coefficient and exact kinematic factor attached to one Fourier graph."""
struct FourierContribution{C<:Number}
    coefficient::C
    kinematic::MomentumPolynomial{ComplexRationals}
end

function Base.isequal(
    a::FourierContribution{C}, b::FourierContribution{C}
) where {C<:Number}
    return isequal(a.coefficient, b.coefficient) && isequal(a.kinematic, b.kinematic)
end
Base.:(==)(a::FourierContribution, b::FourierContribution) = isequal(a, b)
function Base.hash(contribution::FourierContribution, h::UInt)
    return hash(
        FourierContribution, hash(contribution.kinematic, hash(contribution.coefficient, h))
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
    diagrams::Dict{FourierDiagram{S,E1,E2,Nothing},Vector{FourierContribution{C}}}
end

function FourierDiagrams{C,S,E1,E2}() where {C<:Number,S<:Statistics,E1,E2}
    K = FourierDiagram{S,E1,E2,Nothing}
    return FourierDiagrams{C,S,E1,E2}(Dict{K,Vector{FourierContribution{C}}}())
end

function SmallCollections.default(::Type{FourierDiagrams{C,S,E1,E2}}) where {C,S,E1,E2}
    return FourierDiagrams{C,S,E1,E2}()
end

Base.length(collection::FourierDiagrams) = length(collection.diagrams)
Base.isempty(collection::FourierDiagrams) = isempty(collection.diagrams)
Base.iszero(collection::FourierDiagrams) = isempty(collection.diagrams)
Base.iterate(collection::FourierDiagrams) = iterate(collection.diagrams)
Base.iterate(collection::FourierDiagrams, state) = iterate(collection.diagrams, state)
function Base.eltype(::Type{FourierDiagrams{C,S,E1,E2}}) where {C,S,E1,E2}
    return Pair{FourierDiagram{S,E1,E2,Nothing},Vector{FourierContribution{C}}}
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
    coefficient::D,
    kinematic::MomentumPolynomial{ComplexRationals},
) where {C<:Number,D<:Number,S<:Statistics,E1,E2}
    value = _simplify(convert(C, coefficient))
    iszero(value) && return collection

    contributions = get!(collection.diagrams, graph) do
        return FourierContribution{C}[]
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

@inline function _without_derivatives(field::Field{S})::Field{S} where {S<:Statistics}
    return Field{S}(
        field.family,
        field.orientation,
        field.keldysh,
        field.position,
        field.regularisation,
        DerivativeMultiIndex(),
    )
end

function _without_derivatives(edge::Edge{S})::Edge{S} where {S<:Statistics}
    return Edge{S}(
        _without_derivatives(edge.out),
        _without_derivatives(edge.in),
        edge.edgetype,
        Momenta(),
    )
end

@inline function _field_fourier_key(field::Field{S}) where {S<:Statistics}
    return (
        field.family.name,
        slots(field.family.indices),
        Int(field.orientation),
        Int(field.keldysh),
        Int(field.regularisation),
    )
end

@inline function _edge_fourier_key(edge::Edge{S}) where {S<:Statistics}
    return (
        sort_by_position_and_type(edge),
        _field_fourier_key(edge.out),
        _field_fourier_key(edge.in),
        Int(edge.edgetype),
    )
end

@inline function _derivative_tiebreak(edge::Edge{S}) where {S<:Statistics}
    return (edge.out.derivative.orders, edge.in.derivative.orders)
end

@inline function _fourier_source_isless(
    stripped_a::Edge{S}, decorated_a::Edge{S}, stripped_b::Edge{S}, decorated_b::Edge{S}
) where {S<:Statistics}
    key_a = _edge_fourier_key(stripped_a)
    key_b = _edge_fourier_key(stripped_b)
    key_a == key_b || return isless(key_a, key_b)
    derivative_a = _derivative_tiebreak(decorated_a)
    derivative_b = _derivative_tiebreak(decorated_b)
    derivative_a == derivative_b && return false
    return isless(derivative_a, derivative_b)
end

@inline function _relabel_fourier_field(
    field::Field{S}, mapping::Dict{Position,Position}
)::Field{S} where {S<:Statistics}
    p = field.position
    if is_bulk(p) && haskey(mapping, p)
        return Field{S}(
            field.family,
            field.orientation,
            field.keldysh,
            mapping[p],
            field.regularisation,
            field.derivative,
        )
    end
    return field
end

@inline function _relabel_fourier_edge(
    edge::Edge{S}, mapping::Dict{Position,Position}
)::Edge{S} where {S<:Statistics}
    return Edge{S}(
        _relabel_fourier_field(edge.out, mapping),
        _relabel_fourier_field(edge.in, mapping),
        edge.edgetype,
        edge.momenta,
    )
end

function _sort_fourier_source!(
    stripped::Vector{Edge{S}}, decorated::Vector{Edge{S}}
)::Nothing where {S<:Statistics}
    @inbounds for i in 2:length(stripped)
        j = i
        while j > 1 && _fourier_source_isless(
            stripped[j], decorated[j], stripped[j - 1], decorated[j - 1]
        )
            stripped[j - 1], stripped[j] = stripped[j], stripped[j - 1]
            decorated[j - 1], decorated[j] = decorated[j], decorated[j - 1]
            j -= 1
        end
    end
    return nothing
end

@inline function _routing_only_fourier_diagram(
    diagram::Diagram{S,E1,E2}
)::FourierDiagram{S,E1,E2,Nothing} where {S<:Statistics,E1,E2}
    basis, momenta, external_count, loop_count = _diagram_affine_routing(diagram)
    return FourierDiagram{S,E1,E2,Nothing}(
        diagram, basis, momenta, external_count, loop_count, nothing
    )
end

function _canonical_fourier_source(
    diagram::Diagram{S,E1,E2}
)::Tuple{
    FourierDiagram{S,E1,E2,Nothing},
    MomentumPolynomial{ComplexRationals},
    FourierDiagram{S,E1,E2,MomentumPolynomial{ComplexRationals}},
} where {S<:Statistics,E1,E2}
    source_edges = contractions(diagram)
    decorated = Vector{Edge{S}}(undef, E1)
    stripped = Vector{Edge{S}}(undef, E1)
    @inbounds for i in 1:E1
        edge = source_edges[i]
        copied = Edge{S}(edge.out, edge.in, edge.edgetype, Momenta())
        decorated[i] = copied
        stripped[i] = _without_derivatives(copied)
    end

    graph_positions = canonicalization_positions(stripped)
    physical_permutation, _, _, _ = canonicalization_permutations(stripped, graph_positions)
    mapping = make_permutation_dict(physical_permutation, graph_positions, stripped)

    @inbounds for i in 1:E1
        decorated[i] = _relabel_fourier_edge(decorated[i], mapping)
        stripped[i] = _relabel_fourier_edge(stripped[i], mapping)
    end
    _sort_fourier_source!(stripped, decorated)

    stripped_fixed = FixedVector{E1,Edge{S}}(stripped)
    stripped_diagram = Diagram(stripped_fixed, Val(E2))::Diagram{S,E1,E2}
    decorated_fixed = FixedVector{E1,Edge{S}}(decorated)
    decorated_diagram = Diagram{S,E1,E2}(decorated_fixed, topology(stripped_diagram))

    routed = _routing_only_fourier_diagram(stripped_diagram)
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
    diagrams::Diagrams{C,S,E1,E2}
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
    a::FourierDressedPropagator{C,S,O,E1,E2}, b::FourierDressedPropagator{C,S,O,E1,E2}
) where {C,S,O,E1,E2}
    return isequal(a.keldysh, b.keldysh) &&
           isequal(a.retarded, b.retarded) &&
           isequal(a.advanced, b.advanced) &&
           isequal(a.parameter, b.parameter)
end
Base.:(==)(a::FourierDressedPropagator, b::FourierDressedPropagator) = isequal(a, b)
function Base.hash(G::FourierDressedPropagator, h::UInt)
    return hash((G.keldysh, G.retarded, G.advanced, G.parameter), h)
end

function fourier_transform(
    G::DressedPropagator{C,S,O,E1,E2}
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
    a::FourierSelfEnergy{C,S,O,E1,E2}, b::FourierSelfEnergy{C,S,O,E1,E2}
) where {C,S,O,E1,E2}
    return isequal(a.keldysh, b.keldysh) &&
           isequal(a.retarded, b.retarded) &&
           isequal(a.advanced, b.advanced) &&
           isequal(a.parameter, b.parameter)
end
Base.:(==)(a::FourierSelfEnergy, b::FourierSelfEnergy) = isequal(a, b)
function Base.hash(Σ::FourierSelfEnergy, h::UInt)
    return hash((Σ.keldysh, Σ.retarded, Σ.advanced, Σ.parameter), h)
end

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
    length(edges) == SE ||
        error("Fourier self-energy amputation found the wrong bulk-edge count")

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
    G::FourierDressedPropagator{C,S,O,E1,E2}
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

"""
Reject direct transformation of an amputated coordinate-space self-energy.

`SelfEnergy` does not retain provenance proving that external derivative metadata was absent
before amputation. Transform the `DressedPropagator` first, then construct
`SelfEnergy(fourier_transform(G))` so all external momentum factors are preserved exactly.
"""
function fourier_transform(::SelfEnergy)
    return throw(
        ArgumentError(
            "cannot safely Fourier-transform an amputated coordinate SelfEnergy; transform the DressedPropagator first and construct SelfEnergy(fourier_transform(G))",
        ),
    )
end

function _fourier_zero(::Type{FourierDiagrams{C,S,E1,E2}}) where {C,S,E1,E2}
    return FourierDiagrams{C,S,E1,E2}()
end

function matrix(G::FourierDressedPropagator{C,Boson,O,E1,E2}) where {C<:Number,O,E1,E2}
    D = FourierDiagrams{C,Boson,E1,E2}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = G.keldysh
    result[1, 2] = G.retarded
    result[2, 1] = G.advanced
    result[2, 2] = _fourier_zero(D)
    return result
end

function matrix(G::FourierDressedPropagator{C,Fermion,O,E1,E2}) where {C<:Number,O,E1,E2}
    D = FourierDiagrams{C,Fermion,E1,E2}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = G.retarded
    result[1, 2] = G.keldysh
    result[2, 1] = _fourier_zero(D)
    result[2, 2] = G.advanced
    return result
end

function matrix(Σ::FourierSelfEnergy{C,Boson,O,E1,E2}) where {C<:Number,O,E1,E2}
    D = FourierDiagrams{C,Boson,E1,E2}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = _fourier_zero(D)
    result[1, 2] = Σ.advanced
    result[2, 1] = Σ.retarded
    result[2, 2] = Σ.keldysh
    return result
end

function matrix(Σ::FourierSelfEnergy{C,Fermion,O,E1,E2}) where {C<:Number,O,E1,E2}
    D = FourierDiagrams{C,Fermion,E1,E2}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = Σ.retarded
    result[1, 2] = Σ.keldysh
    result[2, 1] = _fourier_zero(D)
    result[2, 2] = Σ.advanced
    return result
end
