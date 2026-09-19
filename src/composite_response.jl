"""
Opaque homogeneous regular HS response for the physical two-body-loss closure.

`polarization` is the χ-target Fourier self-energy generated from the same 2PI skeleton. The
physical pair parameters are named explicitly because the exact response is non-perturbative in
`g` and `γ` and therefore is not an ordinary `ParameterMonomial`. This object represents
`Dreg = D - D0`, with `D⁻¹ = D0⁻¹ - Ω`, without expanding the rational response into a graph
series.
"""
struct TwoBodyLossHSResponse{C<:Number,O,E1,E2}
    polarization::FourierSelfEnergy{C,Boson,O,E1,E2}
    family::FieldFamily{Boson}
    coherent_parameter::Symbol
    loss_parameter::Symbol
end

function TwoBodyLossHSResponse(
    polarization::FourierSelfEnergy{C,Boson,O,E1,E2},
    family::FieldFamily{Boson};
    coherent_parameter::Symbol=:g,
    loss_parameter::Symbol=:γ,
) where {C<:Number,O,E1,E2}
    isequal(target_family(polarization), family) || throw(
        ArgumentError("HS response family must equal the polarization target family")
    )
    coherent_parameter == loss_parameter && throw(
        ArgumentError("coherent and loss parameters must be distinct")
    )
    return TwoBodyLossHSResponse{C,O,E1,E2}(
        polarization, family, coherent_parameter, loss_parameter
    )
end

response_polarization(response::TwoBodyLossHSResponse) = response.polarization
response_family(response::TwoBodyLossHSResponse) = response.family
response_coherent_parameter(response::TwoBodyLossHSResponse) = response.coherent_parameter
response_loss_parameter(response::TwoBodyLossHSResponse) = response.loss_parameter

"""
One canonically routed Fourier graph carrying an opaque composite response.

The edge selected by `response_edge` remains in the coordinate graph only as a routing and
endpoint-provenance placeholder. It must never be interpreted as an independent statistical χ
propagator. Its momentum and Keldysh component specify which component of the opaque regular
response multiplies the remaining physical ψ graph.
"""
struct CompositeFourierDiagram{E1,E2}
    fourier::FourierDiagram{Boson,E1,E2,Nothing}
    response_edge::Int16
end

function CompositeFourierDiagram(
    fourier::FourierDiagram{Boson,E1,E2,Nothing}, response_edge::Integer
) where {E1,E2}
    1 <= response_edge <= E1 || throw(ArgumentError("invalid composite response edge"))
    return CompositeFourierDiagram{E1,E2}(fourier, convert(Int16, response_edge))
end

fourier_diagram(diagram::CompositeFourierDiagram) = diagram.fourier
coordinate_diagram(diagram::CompositeFourierDiagram) = coordinate_diagram(diagram.fourier)
momentum_basis(diagram::CompositeFourierDiagram) = momentum_basis(diagram.fourier)
edge_momenta(diagram::CompositeFourierDiagram) = edge_momenta(diagram.fourier)
external_momentum_count(diagram::CompositeFourierDiagram) =
    external_momentum_count(diagram.fourier)
loop_momentum_count(diagram::CompositeFourierDiagram) = loop_momentum_count(diagram.fourier)
response_edge_index(diagram::CompositeFourierDiagram) = Int(diagram.response_edge)

function _composite_response_edge(diagram::CompositeFourierDiagram)
    return contractions(coordinate_diagram(diagram))[response_edge_index(diagram)]
end

response_family(diagram::CompositeFourierDiagram) = field_family(_composite_response_edge(diagram).out)
response_component(diagram::CompositeFourierDiagram) =
    propagator_type(_composite_response_edge(diagram))
response_momentum(diagram::CompositeFourierDiagram) =
    edge_momenta(diagram)[response_edge_index(diagram)]

function physical_contractions(diagram::CompositeFourierDiagram)
    index = response_edge_index(diagram)
    return Contraction{Boson}[
        edge for (i, edge) in enumerate(contractions(coordinate_diagram(diagram))) if i != index
    ]
end

function Base.isequal(
    a::CompositeFourierDiagram{E1,E2}, b::CompositeFourierDiagram{E1,E2}
) where {E1,E2}
    return isequal(a.fourier, b.fourier) && a.response_edge == b.response_edge
end
Base.:(==)(a::CompositeFourierDiagram, b::CompositeFourierDiagram) = isequal(a, b)
function Base.hash(diagram::CompositeFourierDiagram, h::UInt)
    return hash(
        CompositeFourierDiagram, hash(diagram.response_edge, hash(diagram.fourier, h))
    )
end

"""Collection of routed outer graphs with one opaque response edge each."""
struct CompositeFourierDiagrams{C<:Number,E1,E2}
    diagrams::Dict{CompositeFourierDiagram{E1,E2},Vector{FourierContribution{C}}}
end

function CompositeFourierDiagrams{C,E1,E2}() where {C<:Number,E1,E2}
    K = CompositeFourierDiagram{E1,E2}
    return CompositeFourierDiagrams{C,E1,E2}(
        Dict{K,Vector{FourierContribution{C}}}()
    )
end

Base.length(collection::CompositeFourierDiagrams) = length(collection.diagrams)
Base.isempty(collection::CompositeFourierDiagrams) = isempty(collection.diagrams)
Base.iszero(collection::CompositeFourierDiagrams) = isempty(collection.diagrams)
Base.iterate(collection::CompositeFourierDiagrams) = iterate(collection.diagrams)
Base.iterate(collection::CompositeFourierDiagrams, state) =
    iterate(collection.diagrams, state)
function Base.isequal(
    a::CompositeFourierDiagrams{C,E1,E2}, b::CompositeFourierDiagrams{C,E1,E2}
) where {C<:Number,E1,E2}
    return isequal(a.diagrams, b.diagrams)
end
Base.:(==)(a::CompositeFourierDiagrams, b::CompositeFourierDiagrams) = isequal(a, b)
Base.hash(collection::CompositeFourierDiagrams, h::UInt) = hash(collection.diagrams, h)

function _push_composite_fourier!(
    collection::CompositeFourierDiagrams{C,E1,E2},
    graph::CompositeFourierDiagram{E1,E2},
    coefficient::D,
    kinematic::MomentumPolynomial{ComplexRationals},
) where {C<:Number,D<:Number,E1,E2}
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

"""
Fourier-side self-energy whose regular HS line is an opaque composite response.

Unlike `FourierSelfEnergy`, this object is deliberately not a valid input to the ordinary Wigner
compiler. The response edge still carries χ routing provenance, and a dedicated physical lowering
must eliminate that placeholder into ψ response content before statistical atoms are constructed.
"""
struct CompositeFourierSelfEnergy{C<:Number,O,E1,E2,R}
    keldysh::CompositeFourierDiagrams{C,E1,E2}
    retarded::CompositeFourierDiagrams{C,E1,E2}
    advanced::CompositeFourierDiagrams{C,E1,E2}
    response::R
    target::FieldFamily{Boson}
end

statistics(::CompositeFourierSelfEnergy) = Boson
target_family(self_energy::CompositeFourierSelfEnergy) = self_energy.target
response_family(self_energy::CompositeFourierSelfEnergy) = response_family(self_energy.response)
response_polarization(self_energy::CompositeFourierSelfEnergy) =
    response_polarization(self_energy.response)
response_coherent_parameter(self_energy::CompositeFourierSelfEnergy) =
    response_coherent_parameter(self_energy.response)
response_loss_parameter(self_energy::CompositeFourierSelfEnergy) =
    response_loss_parameter(self_energy.response)
keldysh_component(self_energy::CompositeFourierSelfEnergy) = self_energy.keldysh
retarded_component(self_energy::CompositeFourierSelfEnergy) = self_energy.retarded
advanced_component(self_energy::CompositeFourierSelfEnergy) = self_energy.advanced

function _composite_response_edge_index(
    graph::FourierDiagram{Boson,E1,E2,Nothing}, family::FieldFamily{Boson}
) where {E1,E2}
    found = 0
    for (index, edge) in enumerate(contractions(coordinate_diagram(graph)))
        isequal(field_family(edge.out), family) || continue
        found == 0 || throw(ArgumentError("composite graph contains multiple response edges"))
        found = index
    end
    found == 0 && throw(ArgumentError("composite graph is missing its response edge"))
    return found
end

"""
    twopi_composite_fourier_self_energy(Γ₂, target, response_family; ...)

Construct the production Fourier-side carrier for the physical regular HS response. The target
line is cut while its endpoints are known, physical external ψ segments are attached, and the
outer ψ+response graph is canonically Fourier-routed before amputation. The response-family edge
is retained only as an opaque routing placeholder whose exact `Dreg[Ω[G]]` semantics live in the
attached `TwoBodyLossHSResponse`.

No response denominator is expanded into diagrams and χ is not lowered as a statistical species.
"""
function twopi_composite_fourier_self_energy(
    Γ::TwoPIEffectiveAction{C,Boson,O,E,E2},
    target::FieldFamily{Boson},
    response_family_field::FieldFamily{Boson};
    coherent_parameter::Symbol=:g,
    loss_parameter::Symbol=:γ,
) where {C<:Number,O,E,E2}
    families = field_families(Γ)
    target in families || throw(ArgumentError("target field family is absent from Γ₂"))
    response_family_field in families ||
        throw(ArgumentError("response field family is absent from Γ₂"))
    isequal(target, response_family_field) &&
        throw(ArgumentError("target and response field families must be distinct"))

    polarization = twopi_fourier_self_energy(Γ, response_family_field)
    response = TwoBodyLossHSResponse(
        polarization,
        response_family_field;
        coherent_parameter=coherent_parameter,
        loss_parameter=loss_parameter,
    )

    SE = E - 1
    ST = max_edges(O)
    D = CompositeFourierDiagrams{C,SE,ST}
    components = SmallCollections.SmallDict{3,PropagatorType.T,D}((
        PropagatorType.Advanced => D(),
        PropagatorType.Retarded => D(),
        PropagatorType.Keldysh => D(),
    ))
    derivative_factor = _twopi_derivative_factor(Boson, C)

    for (vacuum, coefficient) in Γ
        contractions_source = twopi_contractions(vacuum)
        for cut_index in eachindex(contractions_source)
            cut = contractions_source[cut_index]
            isequal(field_family(cut.out), target) || continue

            component = _twopi_physical_component(keldysh_index(cut.in), keldysh_index(cut.out))
            iszero(component) && continue

            internal = Contraction{Boson}[
                contractions_source[index] for
                index in eachindex(contractions_source) if index != cut_index
            ]
            _twopi_internal_is_physical(internal) || continue

            source = _twopi_fourier_anchored_diagram(
                contractions_source, cut_index, target, component, Val(ST)
            )
            routed, kinematic, _ = _canonical_fourier_source(source)
            amputated = _amputate_fourier_graph(routed, Val(SE), Val(ST))
            response_index = _composite_response_edge_index(amputated, response_family_field)
            graph = CompositeFourierDiagram(amputated, response_index)

            all(
                edge -> isequal(field_family(edge.out), target), physical_contractions(graph)
            ) || throw(ArgumentError("composite outer graph contains a nonphysical field family"))

            destination = if component == 0x01
                components[PropagatorType.Retarded]
            elseif component == 0x02
                components[PropagatorType.Advanced]
            else
                components[PropagatorType.Keldysh]
            end
            _push_composite_fourier!(
                destination, graph, derivative_factor * coefficient, kinematic
            )
        end
    end

    return CompositeFourierSelfEnergy{C,O,SE,ST,typeof(response)}(
        components[PropagatorType.Keldysh],
        components[PropagatorType.Retarded],
        components[PropagatorType.Advanced],
        response,
        target,
    )
end

function wigner_transform(::CompositeFourierSelfEnergy; gradient_order::Val=Val(0))
    return throw(
        ArgumentError(
            "composite HS response must be eliminated into physical ψ content before ordinary Wigner/statistical lowering",
        ),
    )
end
