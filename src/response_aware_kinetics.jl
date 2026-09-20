"""Exact regular two-body-loss HS response components in KC's symmetric Keldysh basis."""
function _two_body_loss_regular_response(
    g::Number,
    γ::Number,
    ΩK::Number,
    ΩR::Number,
    ΩA::Number,
)
    λ = g - im * γ
    λbar = g + im * γ
    λ2 = g * g + γ * γ
    QR = one(λ) + im * λ * ΩR
    QA = one(λbar) + im * λbar * ΩA

    retarded = -(λ * λ * ΩR) / QR
    advanced = -(λbar * λbar * ΩA) / QA
    keldysh = (
        -λ2 * ΩK +
        2im * γ * (λ * ΩR + λbar * ΩA) -
        2 * γ * λ2 * ΩR * ΩA
    ) / (QR * QA)

    return (; keldysh, retarded, advanced)
end

"""Linear-in-polarization expansion of the exact regular HS response."""
function _two_body_loss_regular_response_linear(
    g::Number,
    γ::Number,
    ΩK::Number,
    ΩR::Number,
    ΩA::Number,
)
    λ = g - im * γ
    λbar = g + im * γ
    λ2 = g * g + γ * γ
    return (
        keldysh=-λ2 * ΩK + 2im * γ * (λ * ΩR + λbar * ΩA),
        retarded=-(λ * λ) * ΩR,
        advanced=-(λbar * λbar) * ΩA,
    )
end

"""Coupling-sector coefficients of the weak regular response."""
function _two_body_loss_regular_response_linear_sectors(
    ΩK::Number, ΩR::Number, ΩA::Number
)
    return (
        g2=(keldysh=-ΩK, retarded=-ΩR, advanced=-ΩA),
        gγ=(
            keldysh=2im * (ΩR + ΩA),
            retarded=2im * ΩR,
            advanced=-2im * ΩA,
        ),
        γ2=(
            keldysh=-ΩK + 2 * ΩR - 2 * ΩA,
            retarded=ΩR,
            advanced=ΩA,
        ),
    )
end

"""
One homogeneous Wigner graph carrying a routed opaque pair response.

The embedded `CompositeFourierDiagram` retains the response edge only as routing provenance.
`physical_contractions` excludes that edge, so downstream physical-line lowering cannot interpret
χ as an independent quasiparticle.
"""
struct ResponseAwareWignerDiagram{E1,E2,Ctx<:AbstractWignerContext}
    composite::CompositeFourierDiagram{E1,E2}
    external_momentum::MomentumVariable
    context::Ctx
end

statistics(::ResponseAwareWignerDiagram) = Boson
gradient_order(::ResponseAwareWignerDiagram) = Val(0)
fourier_diagram(diagram::ResponseAwareWignerDiagram) = fourier_diagram(diagram.composite)
coordinate_diagram(diagram::ResponseAwareWignerDiagram) = coordinate_diagram(diagram.composite)
momentum_basis(diagram::ResponseAwareWignerDiagram) = momentum_basis(diagram.composite)
edge_momenta(diagram::ResponseAwareWignerDiagram) = edge_momenta(diagram.composite)
external_wigner_momentum(diagram::ResponseAwareWignerDiagram) = diagram.external_momentum
wigner_context(diagram::ResponseAwareWignerDiagram) = diagram.context
response_family(diagram::ResponseAwareWignerDiagram) = response_family(diagram.composite)
response_component(diagram::ResponseAwareWignerDiagram) = response_component(diagram.composite)
response_momentum(diagram::ResponseAwareWignerDiagram) = response_momentum(diagram.composite)
physical_contractions(diagram::ResponseAwareWignerDiagram) = physical_contractions(diagram.composite)

function Base.isequal(
    a::ResponseAwareWignerDiagram{E1,E2,Ctx},
    b::ResponseAwareWignerDiagram{E1,E2,Ctx},
) where {E1,E2,Ctx}
    return isequal(a.composite, b.composite) &&
           isequal(a.external_momentum, b.external_momentum) &&
           isequal(a.context, b.context)
end
Base.:(==)(a::ResponseAwareWignerDiagram, b::ResponseAwareWignerDiagram) = isequal(a, b)
function Base.hash(diagram::ResponseAwareWignerDiagram, h::UInt)
    return hash(
        ResponseAwareWignerDiagram,
        hash(diagram.context, hash(diagram.external_momentum, hash(diagram.composite, h))),
    )
end

function _homogeneous_response_wigner_diagram(diagram::CompositeFourierDiagram{E1,E2}) where {E1,E2}
    external_momentum_count(diagram) == 1 || throw(
        ArgumentError("response-aware homogeneous Wigner carrier requires one external momentum"),
    )
    basis = momentum_basis(diagram)
    isempty(basis.variables) && error("response-aware Wigner input has no momentum basis")
    return ResponseAwareWignerDiagram{E1,E2,HomogeneousWignerContext}(
        diagram, basis[1], HomogeneousWignerContext()
    )
end

"""Concrete homogeneous Wigner collection with one opaque response factor per graph."""
struct ResponseAwareWignerDiagrams{C<:Number,E1,E2,Ctx<:AbstractWignerContext}
    diagrams::Dict{ResponseAwareWignerDiagram{E1,E2,Ctx},Vector{WignerContribution{C}}}
    context::Ctx
end

Base.length(collection::ResponseAwareWignerDiagrams) = length(collection.diagrams)
Base.isempty(collection::ResponseAwareWignerDiagrams) = isempty(collection.diagrams)
Base.iszero(collection::ResponseAwareWignerDiagrams) = isempty(collection.diagrams)
Base.iterate(collection::ResponseAwareWignerDiagrams) = iterate(collection.diagrams)
Base.iterate(collection::ResponseAwareWignerDiagrams, state) = iterate(collection.diagrams, state)
wigner_context(collection::ResponseAwareWignerDiagrams) = collection.context
gradient_order(::ResponseAwareWignerDiagrams) = Val(0)

function _homogeneous_response_wigner_diagrams(
    diagrams::CompositeFourierDiagrams{C,E1,E2},
) where {C<:Number,E1,E2}
    context = HomogeneousWignerContext()
    K = ResponseAwareWignerDiagram{E1,E2,HomogeneousWignerContext}
    out = ResponseAwareWignerDiagrams{C,E1,E2,HomogeneousWignerContext}(
        Dict{K,Vector{WignerContribution{C}}}(), context
    )
    for (graph, contributions) in diagrams
        wigner_graph = _homogeneous_response_wigner_diagram(graph)
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

"""
Homogeneous Wigner self-energy with a nested exact HS response.

The response polarization remains a separate pure-ψ object. No χ `KineticLine`, statistical
atom, occupation, or dispersion object is constructed at this boundary.
"""
struct ResponseAwareWignerSelfEnergy{C<:Number,O,E1,E2,R,Ctx<:AbstractWignerContext}
    keldysh::ResponseAwareWignerDiagrams{C,E1,E2,Ctx}
    retarded::ResponseAwareWignerDiagrams{C,E1,E2,Ctx}
    advanced::ResponseAwareWignerDiagrams{C,E1,E2,Ctx}
    response::R
    target::FieldFamily{Boson}
    context::Ctx
end

statistics(::ResponseAwareWignerSelfEnergy) = Boson
order(::ResponseAwareWignerSelfEnergy{C,O}) where {C,O} = O
target_family(self_energy::ResponseAwareWignerSelfEnergy) = self_energy.target
response_family(self_energy::ResponseAwareWignerSelfEnergy) = response_family(self_energy.response)
response_polarization(self_energy::ResponseAwareWignerSelfEnergy) = response_polarization(self_energy.response)
response_coherent_parameter(self_energy::ResponseAwareWignerSelfEnergy) = response_coherent_parameter(self_energy.response)
response_loss_parameter(self_energy::ResponseAwareWignerSelfEnergy) = response_loss_parameter(self_energy.response)
keldysh_component(self_energy::ResponseAwareWignerSelfEnergy) = self_energy.keldysh
retarded_component(self_energy::ResponseAwareWignerSelfEnergy) = self_energy.retarded
advanced_component(self_energy::ResponseAwareWignerSelfEnergy) = self_energy.advanced
gradient_order(::ResponseAwareWignerSelfEnergy) = Val(0)
wigner_context(self_energy::ResponseAwareWignerSelfEnergy) = self_energy.context

function _response_aware_wigner_transform(
    self_energy::CompositeFourierSelfEnergy{C,O,E1,E2,R}, ::Val{0}
) where {C<:Number,O,E1,E2,R}
    context = HomogeneousWignerContext()
    return ResponseAwareWignerSelfEnergy{C,O,E1,E2,R,HomogeneousWignerContext}(
        _homogeneous_response_wigner_diagrams(self_energy.keldysh),
        _homogeneous_response_wigner_diagrams(self_energy.retarded),
        _homogeneous_response_wigner_diagrams(self_energy.advanced),
        self_energy.response,
        target_family(self_energy),
        context,
    )
end

function _response_aware_wigner_transform(
    ::CompositeFourierSelfEnergy, gradient_order::Val{G}
) where {G}
    return _unsupported_wigner_gradient(gradient_order)
end

function response_aware_wigner_transform(
    self_energy::CompositeFourierSelfEnergy; gradient_order::Val=Val(0)
)
    return _response_aware_wigner_transform(self_energy, gradient_order)
end
