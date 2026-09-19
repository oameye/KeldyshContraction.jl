"""
Exact KC-native regular HS response components for the physical two-body-loss closure.

The inputs `ΩK`, `ΩR`, and `ΩA` are the KC-native polarization components at one fixed
frequency-momentum argument. This helper evaluates the closed Dyson response without expanding
its denominators. It is intentionally specialized to the physical pair coupling `λ = g - im*γ`
rather than introducing a generic symbolic rational algebra.
"""
function _two_body_loss_regular_response_components(
    g::Real, γ::Real, ΩK::Number, ΩR::Number, ΩA::Number
)
    λ = complex(g, -γ)
    λbar = complex(g, γ)
    λ2 = λ * λbar
    QR = 1 + im * λ * ΩR
    QA = 1 + im * λbar * ΩA
    denominator = QR * QA

    retarded = -(λ^2) * ΩR / QR
    advanced = -(λbar^2) * ΩA / QA
    keldysh =
        (
            -λ2 * ΩK +
            2im * γ * (λ * ΩR + λbar * ΩA) -
            2 * γ * λ2 * ΩR * ΩA
        ) / denominator
    full_keldysh = (-2 * γ - λ2 * ΩK) / denominator

    return (; keldysh, retarded, advanced, full_keldysh)
end

"""
Response after homogeneous Wigner transformation, with the auxiliary χ target removed.

The three polarization components contain only physical ψ propagator lines. `parameter` retains
the formal 2PI skeleton bookkeeping monomial; the physical non-perturbative pair parameters are
stored separately as names because the resummed response is not a `ParameterMonomial`.
"""
struct TwoBodyLossWignerHSResponse{
    C<:Number,O,E1,E2,G,Ctx<:AbstractWignerContext
}
    keldysh::WignerDiagrams{C,Boson,E1,E2,G,Ctx}
    retarded::WignerDiagrams{C,Boson,E1,E2,G,Ctx}
    advanced::WignerDiagrams{C,Boson,E1,E2,G,Ctx}
    parameter::ParameterMonomial
    physical_family::FieldFamily{Boson}
    coherent_parameter::Symbol
    loss_parameter::Symbol
    context::Ctx
end

order(::TwoBodyLossWignerHSResponse{C,O}) where {C,O} = O
statistics(::TwoBodyLossWignerHSResponse) = Boson
parameters(response::TwoBodyLossWignerHSResponse) = response.parameter
wigner_context(response::TwoBodyLossWignerHSResponse) = response.context
gradient_order(::TwoBodyLossWignerHSResponse{C,O,E1,E2,G}) where {C,O,E1,E2,G} = Val(G)
keldysh_component(response::TwoBodyLossWignerHSResponse) = response.keldysh
retarded_component(response::TwoBodyLossWignerHSResponse) = response.retarded
advanced_component(response::TwoBodyLossWignerHSResponse) = response.advanced
response_physical_family(response::TwoBodyLossWignerHSResponse) = response.physical_family
response_coherent_parameter(response::TwoBodyLossWignerHSResponse) = response.coherent_parameter
response_loss_parameter(response::TwoBodyLossWignerHSResponse) = response.loss_parameter

function _assert_response_wigner_family(
    diagrams::WignerDiagrams{C,Boson,E1,E2,G,Ctx}, family::FieldFamily{Boson}
) where {C<:Number,E1,E2,G,Ctx<:AbstractWignerContext}
    for (graph, _) in diagrams
        for edge in contractions(coordinate_diagram(graph))
            isequal(field_family(edge.out), family) || throw(
                ArgumentError("response polarization contains a nonphysical field family"),
            )
        end
    end
    return nothing
end

function _wigner_hs_response(
    response::TwoBodyLossHSResponse{C,O,E1,E2}, physical_family::FieldFamily{Boson}
) where {C<:Number,O,E1,E2}
    polarization = wigner_transform(response_polarization(response); gradient_order=Val(0))
    _assert_response_wigner_family(polarization.keldysh, physical_family)
    _assert_response_wigner_family(polarization.retarded, physical_family)
    _assert_response_wigner_family(polarization.advanced, physical_family)

    return TwoBodyLossWignerHSResponse{C,O,E1,E2,0,HomogeneousWignerContext}(
        polarization.keldysh,
        polarization.retarded,
        polarization.advanced,
        polarization.parameter,
        physical_family,
        response_coherent_parameter(response),
        response_loss_parameter(response),
        polarization.context,
    )
end

"""
One homogeneous response-aware Wigner graph.

The χ placeholder from `CompositeFourierDiagram` has been removed completely. The remaining
edges and routed momenta are physical ψ propagators. The exact pair response is represented only
by its K/R/A component and routed pair momentum.
"""
struct CompositeWignerDiagram{E1,E2,G,Ctx<:AbstractWignerContext}
    physical_edges::FixedVector{E1,Edge{Boson}}
    physical_momenta::FixedVector{E1,LinearMomentum}
    topology::FixedVector{E2,Int}
    basis::MomentumBasis
    external_momentum::MomentumVariable
    response_component::PropagatorType.T
    response_momentum::LinearMomentum
    context::Ctx
end

statistics(::CompositeWignerDiagram) = Boson
gradient_order(::CompositeWignerDiagram{E1,E2,G}) where {E1,E2,G} = Val(G)
wigner_context(diagram::CompositeWignerDiagram) = diagram.context
momentum_basis(diagram::CompositeWignerDiagram) = diagram.basis
external_wigner_momentum(diagram::CompositeWignerDiagram) = diagram.external_momentum
physical_contractions(diagram::CompositeWignerDiagram) = diagram.physical_edges
physical_momenta(diagram::CompositeWignerDiagram) = diagram.physical_momenta
response_component(diagram::CompositeWignerDiagram) = diagram.response_component
response_momentum(diagram::CompositeWignerDiagram) = diagram.response_momentum

function Base.isequal(
    a::CompositeWignerDiagram{E1,E2,G,Ctx}, b::CompositeWignerDiagram{E1,E2,G,Ctx}
) where {E1,E2,G,Ctx}
    return isequal(a.physical_edges, b.physical_edges) &&
           isequal(a.physical_momenta, b.physical_momenta) &&
           _kinetic_topology_isequal(a.topology, b.topology) &&
           isequal(a.basis, b.basis) &&
           isequal(a.external_momentum, b.external_momentum) &&
           a.response_component === b.response_component &&
           isequal(a.response_momentum, b.response_momentum) &&
           isequal(a.context, b.context)
end
Base.:(==)(a::CompositeWignerDiagram, b::CompositeWignerDiagram) = isequal(a, b)
function Base.hash(diagram::CompositeWignerDiagram{E1,E2,G,Ctx}, h::UInt) where {E1,E2,G,Ctx}
    h = hash(CompositeWignerDiagram, h)
    h = hash(diagram.physical_edges, h)
    h = hash(diagram.physical_momenta, h)
    h = _kinetic_topology_hash(diagram.topology, h)
    h = hash(diagram.basis, h)
    h = hash(diagram.external_momentum, h)
    h = hash(diagram.response_component, h)
    h = hash(diagram.response_momentum, h)
    return hash(diagram.context, h)
end

"""Collection of response-aware homogeneous Wigner graphs."""
struct CompositeWignerDiagrams{
    C<:Number,E1,E2,G,Ctx<:AbstractWignerContext
}
    diagrams::Dict{CompositeWignerDiagram{E1,E2,G,Ctx},Vector{WignerContribution{C}}}
    context::Ctx
end

function CompositeWignerDiagrams{C,E1,E2,G,Ctx}(
    context::Ctx
) where {C<:Number,E1,E2,G,Ctx<:AbstractWignerContext}
    K = CompositeWignerDiagram{E1,E2,G,Ctx}
    return CompositeWignerDiagrams{C,E1,E2,G,Ctx}(
        Dict{K,Vector{WignerContribution{C}}}(), context
    )
end

Base.length(collection::CompositeWignerDiagrams) = length(collection.diagrams)
Base.isempty(collection::CompositeWignerDiagrams) = isempty(collection.diagrams)
Base.iszero(collection::CompositeWignerDiagrams) = isempty(collection.diagrams)
Base.iterate(collection::CompositeWignerDiagrams) = iterate(collection.diagrams)
function Base.iterate(collection::CompositeWignerDiagrams, state)
    return iterate(collection.diagrams, state)
end
wigner_context(collection::CompositeWignerDiagrams) = collection.context

function _homogeneous_composite_wigner_diagram(
    graph::CompositeFourierDiagram{E1,E2}
) where {E1,E2}
    external_momentum_count(graph) == 1 || throw(
        ArgumentError("response-aware Wigner lowering requires exactly one external momentum"),
    )
    physical_count = E1 - 1
    response_index = response_edge_index(graph)
    coordinate = coordinate_diagram(graph)
    source_edges = contractions(coordinate)
    routed = edge_momenta(graph)

    edges = Edge{Boson}[]
    momenta = LinearMomentum[]
    sizehint!(edges, physical_count)
    sizehint!(momenta, physical_count)
    for i in 1:E1
        i == response_index && continue
        push!(edges, source_edges[i])
        push!(momenta, routed[i])
    end
    length(edges) == physical_count || error("invalid physical edge count in composite graph")

    basis = momentum_basis(graph)
    isempty(basis.variables) && error("response-aware Wigner input has no momentum basis")
    context = HomogeneousWignerContext()
    return CompositeWignerDiagram{physical_count,E2,0,HomogeneousWignerContext}(
        FixedVector{physical_count,Edge{Boson}}(edges),
        FixedVector{physical_count,LinearMomentum}(momenta),
        topology(coordinate),
        basis,
        basis[1],
        response_component(graph),
        response_momentum(graph),
        context,
    )
end

function _homogeneous_composite_wigner_diagrams(
    diagrams::CompositeFourierDiagrams{C,E1,E2}
) where {C<:Number,E1,E2}
    physical_count = E1 - 1
    context = HomogeneousWignerContext()
    out = CompositeWignerDiagrams{
        C,physical_count,E2,0,HomogeneousWignerContext
    }(context)
    for (graph, contributions) in diagrams
        wigner_graph = _homogeneous_composite_wigner_diagram(graph)
        values = Vector{WignerContribution{C}}(undef, length(contributions))
        @inbounds for i in eachindex(contributions)
            contribution = contributions[i]
            values[i] = WignerContribution{C}(
                contribution.coefficient, contribution.kinematic
            )
        end
        out.diagrams[wigner_graph] = values
    end
    return out
end

"""Physical response-aware self-energy after homogeneous Wigner lowering."""
struct CompositeWignerSelfEnergy{
    C<:Number,O,E1,E2,G,Ctx<:AbstractWignerContext,R
}
    keldysh::CompositeWignerDiagrams{C,E1,E2,G,Ctx}
    retarded::CompositeWignerDiagrams{C,E1,E2,G,Ctx}
    advanced::CompositeWignerDiagrams{C,E1,E2,G,Ctx}
    response::R
    target::FieldFamily{Boson}
    context::Ctx
end

order(::CompositeWignerSelfEnergy{C,O}) where {C,O} = O
statistics(::CompositeWignerSelfEnergy) = Boson
target_family(self_energy::CompositeWignerSelfEnergy) = self_energy.target
gradient_order(::CompositeWignerSelfEnergy{C,O,E1,E2,G}) where {C,O,E1,E2,G} = Val(G)
wigner_context(self_energy::CompositeWignerSelfEnergy) = self_energy.context
keldysh_component(self_energy::CompositeWignerSelfEnergy) = self_energy.keldysh
retarded_component(self_energy::CompositeWignerSelfEnergy) = self_energy.retarded
advanced_component(self_energy::CompositeWignerSelfEnergy) = self_energy.advanced
response_polarization(self_energy::CompositeWignerSelfEnergy) = self_energy.response
response_coherent_parameter(self_energy::CompositeWignerSelfEnergy) =
    response_coherent_parameter(self_energy.response)
response_loss_parameter(self_energy::CompositeWignerSelfEnergy) =
    response_loss_parameter(self_energy.response)

"""
    response_wigner_transform(Σ::CompositeFourierSelfEnergy; gradient_order=Val(0))

Dedicated physical Wigner lowering for the C4 composite response carrier. Unlike ordinary
`wigner_transform`, this operation removes the χ routing placeholder and preserves only the
response component/momentum together with physical ψ propagator content.
"""
function response_wigner_transform(
    self_energy::CompositeFourierSelfEnergy{C,O,E1,E2,R}; gradient_order::Val=Val(0)
) where {C<:Number,O,E1,E2,R}
    return _response_wigner_transform(self_energy, gradient_order)
end

function _response_wigner_transform(
    self_energy::CompositeFourierSelfEnergy{C,O,E1,E2,R}, ::Val{0}
) where {C<:Number,O,E1,E2,R}
    physical_count = E1 - 1
    response = _wigner_hs_response(self_energy.response, target_family(self_energy))
    context = HomogeneousWignerContext()
    return CompositeWignerSelfEnergy{
        C,O,physical_count,E2,0,HomogeneousWignerContext,typeof(response)
    }(
        _homogeneous_composite_wigner_diagrams(self_energy.keldysh),
        _homogeneous_composite_wigner_diagrams(self_energy.retarded),
        _homogeneous_composite_wigner_diagrams(self_energy.advanced),
        response,
        target_family(self_energy),
        context,
    )
end

function _response_wigner_transform(::CompositeFourierSelfEnergy, gradient_order::Val{G}) where {G}
    return _unsupported_wigner_gradient(gradient_order)
end

"""Pure-ψ polarization after spectral/statistical lowering."""
struct TwoBodyLossKineticHSResponse{
    C<:Number,O,E1,E2,G,Ctx<:AbstractWignerContext
}
    keldysh::KineticExpression{C,Boson,E1,E2,G,Ctx}
    retarded::KineticExpression{C,Boson,E1,E2,G,Ctx}
    advanced::KineticExpression{C,Boson,E1,E2,G,Ctx}
    parameter::ParameterMonomial
    physical_family::FieldFamily{Boson}
    coherent_parameter::Symbol
    loss_parameter::Symbol
    context::Ctx
end

order(::TwoBodyLossKineticHSResponse{C,O}) where {C,O} = O
statistics(::TwoBodyLossKineticHSResponse) = Boson
parameters(response::TwoBodyLossKineticHSResponse) = response.parameter
wigner_context(response::TwoBodyLossKineticHSResponse) = response.context
gradient_order(::TwoBodyLossKineticHSResponse{C,O,E1,E2,G}) where {C,O,E1,E2,G} = Val(G)
keldysh_component(response::TwoBodyLossKineticHSResponse) = response.keldysh
retarded_component(response::TwoBodyLossKineticHSResponse) = response.retarded
advanced_component(response::TwoBodyLossKineticHSResponse) = response.advanced
response_physical_family(response::TwoBodyLossKineticHSResponse) = response.physical_family
response_coherent_parameter(response::TwoBodyLossKineticHSResponse) = response.coherent_parameter
response_loss_parameter(response::TwoBodyLossKineticHSResponse) = response.loss_parameter

function _assert_response_kinetic_family(
    expression::KineticExpression{C,Boson,E1,E2,G,Ctx}, family::FieldFamily{Boson}
) where {C<:Number,E1,E2,G,Ctx<:AbstractWignerContext}
    for (term, _) in expression
        for line in kinetic_lines(term)
            isequal(line.family, family) || throw(
                ArgumentError("response kinetic polarization contains a nonphysical field family"),
            )
        end
    end
    return nothing
end

function _kinetic_hs_response(
    response::TwoBodyLossWignerHSResponse{C,O,E1,E2,0,Ctx}
) where {C<:Number,O,E1,E2,Ctx<:AbstractWignerContext}
    K = kinetic_expression(response.keldysh)
    R = kinetic_expression(response.retarded)
    A = kinetic_expression(response.advanced)
    family = response.physical_family
    _assert_response_kinetic_family(K, family)
    _assert_response_kinetic_family(R, family)
    _assert_response_kinetic_family(A, family)
    D = kinetic_coefficient_type(C)
    return TwoBodyLossKineticHSResponse{D,O,E1,E2,0,Ctx}(
        K,
        R,
        A,
        response.parameter,
        family,
        response.coherent_parameter,
        response.loss_parameter,
        response.context,
    )
end

"""Opaque exact response factor attached to one physical kinetic term."""
struct TwoBodyLossResponseFactor
    component::PropagatorType.T
    momentum::LinearMomentum
end

response_component(factor::TwoBodyLossResponseFactor) = factor.component
response_momentum(factor::TwoBodyLossResponseFactor) = factor.momentum
Base.isequal(a::TwoBodyLossResponseFactor, b::TwoBodyLossResponseFactor) =
    a.component === b.component && isequal(a.momentum, b.momentum)
Base.:(==)(a::TwoBodyLossResponseFactor, b::TwoBodyLossResponseFactor) = isequal(a, b)
Base.hash(factor::TwoBodyLossResponseFactor, h::UInt) =
    hash(TwoBodyLossResponseFactor, hash(factor.momentum, hash(factor.component, h)))

"""One physical ψ kinetic term multiplied by one exact nested pair response."""
struct CompositeKineticTerm{E1,E2}
    monomial::KineticMonomial{Boson,E1}
    topology::FixedVector{E2,Int}
    basis::MomentumBasis
    external_momentum::MomentumVariable
    response::TwoBodyLossResponseFactor
    kinematic::MomentumPolynomial{ComplexRationals}
end

statistics(::CompositeKineticTerm) = Boson
kinetic_lines(term::CompositeKineticTerm) = kinetic_lines(term.monomial)
kinematic_factor(term::CompositeKineticTerm) = term.kinematic
momentum_basis(term::CompositeKineticTerm) = term.basis
external_wigner_momentum(term::CompositeKineticTerm) = term.external_momentum
response_component(term::CompositeKineticTerm) = response_component(term.response)
response_momentum(term::CompositeKineticTerm) = response_momentum(term.response)

function Base.isequal(
    a::CompositeKineticTerm{E1,E2}, b::CompositeKineticTerm{E1,E2}
) where {E1,E2}
    return isequal(a.monomial, b.monomial) &&
           _kinetic_topology_isequal(a.topology, b.topology) &&
           isequal(a.basis, b.basis) &&
           isequal(a.external_momentum, b.external_momentum) &&
           isequal(a.response, b.response) &&
           isequal(a.kinematic, b.kinematic)
end
Base.:(==)(a::CompositeKineticTerm, b::CompositeKineticTerm) = isequal(a, b)
function Base.hash(term::CompositeKineticTerm{E1,E2}, h::UInt) where {E1,E2}
    h = hash(CompositeKineticTerm, h)
    h = hash(term.monomial, h)
    h = _kinetic_topology_hash(term.topology, h)
    h = hash(term.basis, h)
    h = hash(term.external_momentum, h)
    h = hash(term.response, h)
    return hash(term.kinematic, h)
end

"""Numeric linear combination of physical kinetic terms with exact nested response factors."""
struct CompositeKineticExpression{
    C<:Number,E1,E2,G,Ctx<:AbstractWignerContext
}
    terms::Dict{CompositeKineticTerm{E1,E2},C}
    context::Ctx
end

function CompositeKineticExpression{C,E1,E2,G,Ctx}(
    context::Ctx
) where {C<:Number,E1,E2,G,Ctx<:AbstractWignerContext}
    return CompositeKineticExpression{C,E1,E2,G,Ctx}(
        Dict{CompositeKineticTerm{E1,E2},C}(), context
    )
end

statistics(::CompositeKineticExpression) = Boson
gradient_order(::CompositeKineticExpression{C,E1,E2,G}) where {C,E1,E2,G} = Val(G)
wigner_context(expression::CompositeKineticExpression) = expression.context
Base.length(expression::CompositeKineticExpression) = length(expression.terms)
Base.isempty(expression::CompositeKineticExpression) = isempty(expression.terms)
Base.iszero(expression::CompositeKineticExpression) = isempty(expression.terms)
Base.iterate(expression::CompositeKineticExpression) = iterate(expression.terms)
Base.iterate(expression::CompositeKineticExpression, state) = iterate(expression.terms, state)

function Base.push!(
    expression::CompositeKineticExpression{C,E1,E2,G,Ctx},
    term::CompositeKineticTerm{E1,E2},
    coefficient::Number,
) where {C<:Number,E1,E2,G,Ctx<:AbstractWignerContext}
    value = _simplify(convert(C, coefficient))
    iszero(value) && return expression
    if haskey(expression.terms, term)
        combined = _simplify(expression.terms[term] + value)
        if iszero(combined)
            delete!(expression.terms, term)
        else
            expression.terms[term] = combined
        end
    else
        expression.terms[term] = value
    end
    return expression
end

function Base.:+(
    a::CompositeKineticExpression{C1,E1,E2,G,Ctx},
    b::CompositeKineticExpression{C2,E1,E2,G,Ctx},
) where {C1<:Number,C2<:Number,E1,E2,G,Ctx<:AbstractWignerContext}
    isequal(a.context, b.context) || throw(
        ArgumentError("cannot add composite kinetic expressions with different Wigner contexts"),
    )
    D = promote_type(C1, C2)
    out = CompositeKineticExpression{D,E1,E2,G,Ctx}(a.context)
    for (term, coefficient) in a
        push!(out, term, convert(D, coefficient))
    end
    for (term, coefficient) in b
        push!(out, term, convert(D, coefficient))
    end
    return out
end

function Base.:*(
    prefactor::P, expression::CompositeKineticExpression{C,E1,E2,G,Ctx}
) where {P<:Number,C<:Number,E1,E2,G,Ctx<:AbstractWignerContext}
    D = promote_type(P, C)
    out = CompositeKineticExpression{D,E1,E2,G,Ctx}(expression.context)
    p = convert(D, prefactor)
    for (term, coefficient) in expression
        push!(out, term, p * convert(D, coefficient))
    end
    return out
end
Base.:*(expression::CompositeKineticExpression, prefactor::Number) = prefactor * expression
function Base.:-(expression::CompositeKineticExpression{C}) where {C<:Number}
    return -one(C) * expression
end
Base.:-(a::CompositeKineticExpression, b::CompositeKineticExpression) = a + (-b)

function _lower_composite_kinetic_term(
    graph::CompositeWignerDiagram{E1,E2,0,Ctx},
    contribution::WignerContribution{C},
    ::Type{D},
) where {C<:Number,D<:Number,E1,E2,Ctx<:AbstractWignerContext}
    lines = Vector{KineticLine{Boson}}(undef, E1)
    keldysh_count = Int8(0)
    @inbounds for i in 1:E1
        line, is_statistical =
            _lower_kinetic_line(graph.physical_edges[i], graph.physical_momenta[i])
        lines[i] = line
        keldysh_count += is_statistical
    end

    term = CompositeKineticTerm{E1,E2}(
        KineticMonomial(lines, Val(E1)),
        graph.topology,
        graph.basis,
        graph.external_momentum,
        TwoBodyLossResponseFactor(graph.response_component, graph.response_momentum),
        contribution.kinematic,
    )
    phase = convert(D, (-im)^Int(keldysh_count))
    coefficient = _simplify(convert(D, contribution.coefficient) * phase)
    return term, coefficient
end

function kinetic_expression(
    diagrams::CompositeWignerDiagrams{C,E1,E2,0,Ctx}
) where {C<:Number,E1,E2,Ctx<:AbstractWignerContext}
    D = kinetic_coefficient_type(C)
    out = CompositeKineticExpression{D,E1,E2,0,Ctx}(wigner_context(diagrams))
    for (graph, contributions) in diagrams
        for contribution in contributions
            term, coefficient = _lower_composite_kinetic_term(graph, contribution, D)
            push!(out, term, coefficient)
        end
    end
    return out
end

function kinetic_expression(
    ::CompositeWignerDiagrams{C,E1,E2,G,Ctx}
) where {C<:Number,E1,E2,G,Ctx<:AbstractWignerContext}
    return throw(
        ArgumentError(
            "response-aware kinetic lowering supports only Wigner gradient order zero; got $G",
        ),
    )
end

"""Response-aware physical ψ self-energy with no auxiliary-field kinetic lines."""
struct CompositeKineticSelfEnergy{
    C<:Number,O,E1,E2,G,Ctx<:AbstractWignerContext,R
}
    keldysh::CompositeKineticExpression{C,E1,E2,G,Ctx}
    retarded::CompositeKineticExpression{C,E1,E2,G,Ctx}
    advanced::CompositeKineticExpression{C,E1,E2,G,Ctx}
    response::R
    target::FieldFamily{Boson}
    context::Ctx
end

order(::CompositeKineticSelfEnergy{C,O}) where {C,O} = O
statistics(::CompositeKineticSelfEnergy) = Boson
target_family(self_energy::CompositeKineticSelfEnergy) = self_energy.target
gradient_order(::CompositeKineticSelfEnergy{C,O,E1,E2,G}) where {C,O,E1,E2,G} = Val(G)
wigner_context(self_energy::CompositeKineticSelfEnergy) = self_energy.context
keldysh_component(self_energy::CompositeKineticSelfEnergy) = self_energy.keldysh
retarded_component(self_energy::CompositeKineticSelfEnergy) = self_energy.retarded
advanced_component(self_energy::CompositeKineticSelfEnergy) = self_energy.advanced
response_polarization(self_energy::CompositeKineticSelfEnergy) = self_energy.response
response_coherent_parameter(self_energy::CompositeKineticSelfEnergy) =
    response_coherent_parameter(self_energy.response)
response_loss_parameter(self_energy::CompositeKineticSelfEnergy) =
    response_loss_parameter(self_energy.response)

function kinetic_expression(
    self_energy::CompositeWignerSelfEnergy{C,O,E1,E2,0,Ctx,R}
) where {C<:Number,O,E1,E2,Ctx<:AbstractWignerContext,R}
    response = _kinetic_hs_response(self_energy.response)
    D = kinetic_coefficient_type(C)
    return CompositeKineticSelfEnergy{D,O,E1,E2,0,Ctx,typeof(response)}(
        kinetic_expression(self_energy.keldysh),
        kinetic_expression(self_energy.retarded),
        kinetic_expression(self_energy.advanced),
        response,
        self_energy.target,
        self_energy.context,
    )
end

function kinetic_expression(
    ::CompositeWignerSelfEnergy{C,O,E1,E2,G,Ctx,R}
) where {C<:Number,O,E1,E2,G,Ctx<:AbstractWignerContext,R}
    return throw(
        ArgumentError(
            "response-aware kinetic lowering supports only Wigner gradient order zero; got $G",
        ),
    )
end

"""
Exact off-shell Kadanoff--Baym collision carrier for a resummed pair response.

The collision remains affine in the external physical statistical distribution. Its internal
response factors are opaque `Dreg[Ω[G]]` evaluations and therefore are not sent to the ordinary
quasiparticle frequency-reduction pipeline.
"""
struct CompositeOffShellCollisionExpression{
    C<:Number,O,E1,E2,G,Ctx<:AbstractWignerContext,R
}
    offset::CompositeKineticExpression{C,E1,E2,G,Ctx}
    distribution_coefficient::CompositeKineticExpression{C,E1,E2,G,Ctx}
    response::R
    target::FieldFamily{Boson}
    context::Ctx
end

order(::CompositeOffShellCollisionExpression{C,O}) where {C,O} = O
statistics(::CompositeOffShellCollisionExpression) = Boson
target_family(collision::CompositeOffShellCollisionExpression) = collision.target
gradient_order(::CompositeOffShellCollisionExpression{C,O,E1,E2,G}) where {C,O,E1,E2,G} = Val(G)
wigner_context(collision::CompositeOffShellCollisionExpression) = collision.context
collision_offset(collision::CompositeOffShellCollisionExpression) = collision.offset
function collision_distribution_coefficient(collision::CompositeOffShellCollisionExpression)
    return collision.distribution_coefficient
end
response_polarization(collision::CompositeOffShellCollisionExpression) = collision.response
response_coherent_parameter(collision::CompositeOffShellCollisionExpression) =
    response_coherent_parameter(collision.response)
response_loss_parameter(collision::CompositeOffShellCollisionExpression) =
    response_loss_parameter(collision.response)

function off_shell_collision_expression(
    self_energy::CompositeKineticSelfEnergy{C,O,E1,E2,0,Ctx,R}
) where {C<:Number,O,E1,E2,Ctx<:AbstractWignerContext,R}
    D = kinetic_coefficient_type(C)
    offset = convert(D, im) * self_energy.keldysh
    distribution_coefficient = -convert(D, im) * (self_energy.retarded - self_energy.advanced)
    return CompositeOffShellCollisionExpression{D,O,E1,E2,0,Ctx,R}(
        offset,
        distribution_coefficient,
        self_energy.response,
        self_energy.target,
        self_energy.context,
    )
end

function off_shell_collision_expression(
    ::CompositeKineticSelfEnergy{C,O,E1,E2,G,Ctx,R}
) where {C<:Number,O,E1,E2,G,Ctx<:AbstractWignerContext,R}
    return throw(
        ArgumentError(
            "response-aware off-shell collision construction supports only Wigner gradient order zero; got $G",
        ),
    )
end

function spectral_dispersive_collision(::CompositeOffShellCollisionExpression)
    return throw(
        ArgumentError(
            "resummed pair response has no assumed quasiparticle shell; use a controlled weak-response expansion or a dedicated response spectral backend before ordinary frequency reduction",
        ),
    )
end
