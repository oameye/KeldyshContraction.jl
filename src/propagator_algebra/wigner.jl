abstract type AbstractWignerContext end

"""
Center-coordinate context for the initial homogeneous Wigner representation.

`center_coordinate === nothing` records that no explicit center-coordinate dependence is
represented in the supported gradient-order-zero regime. Future inhomogeneous contexts can
carry center-coordinate metadata without changing the Wigner object families.
"""
struct HomogeneousWignerContext <: AbstractWignerContext
    center_coordinate::Nothing
end

HomogeneousWignerContext() = HomogeneousWignerContext(nothing)

function Base.isequal(a::HomogeneousWignerContext, b::HomogeneousWignerContext)
    return isequal(a.center_coordinate, b.center_coordinate)
end
Base.:(==)(a::HomogeneousWignerContext, b::HomogeneousWignerContext) = isequal(a, b)
function Base.hash(context::HomogeneousWignerContext, h::UInt)
    return hash(HomogeneousWignerContext, hash(context.center_coordinate, h))
end

"""
One Wigner-space diagram.

At gradient order zero the exact Fourier graph is preserved byte-for-byte. The Wigner layer
adds an explicit external relative-coordinate momentum, a center-coordinate context, and the
statically visible gradient order without re-routing or re-deriving derivative kinematics.
"""
struct WignerDiagram{S<:Statistics,E1,E2,G,K,Ctx<:AbstractWignerContext}
    fourier::FourierDiagram{S,E1,E2,K}
    external_momentum::MomentumVariable
    context::Ctx
end

statistics(::WignerDiagram{S}) where {S} = S
gradient_order(::WignerDiagram{S,E1,E2,G}) where {S,E1,E2,G} = Val(G)
fourier_diagram(diagram::WignerDiagram) = diagram.fourier
external_wigner_momentum(diagram::WignerDiagram) = diagram.external_momentum
wigner_context(diagram::WignerDiagram) = diagram.context
coordinate_diagram(diagram::WignerDiagram) = coordinate_diagram(diagram.fourier)
momentum_basis(diagram::WignerDiagram) = momentum_basis(diagram.fourier)
edge_momenta(diagram::WignerDiagram) = edge_momenta(diagram.fourier)
external_momentum_count(diagram::WignerDiagram) = external_momentum_count(diagram.fourier)
loop_momentum_count(diagram::WignerDiagram) = loop_momentum_count(diagram.fourier)
kinematic_factor(diagram::WignerDiagram) = kinematic_factor(diagram.fourier)

function Base.isequal(
    a::WignerDiagram{S,E1,E2,G,K,Ctx}, b::WignerDiagram{S,E1,E2,G,K,Ctx}
) where {S,E1,E2,G,K,Ctx}
    return isequal(a.fourier, b.fourier) &&
           isequal(a.external_momentum, b.external_momentum) &&
           isequal(a.context, b.context)
end
Base.:(==)(a::WignerDiagram, b::WignerDiagram) = isequal(a, b)
function Base.hash(diagram::WignerDiagram, h::UInt)
    h = hash(WignerDiagram, h)
    h = hash(diagram.fourier, h)
    h = hash(diagram.external_momentum, h)
    return hash(diagram.context, h)
end

"""Numeric coefficient and exact kinematic factor attached to one Wigner graph."""
struct WignerContribution{C<:Number}
    coefficient::C
    kinematic::MomentumPolynomial{ComplexRationals}
end

function Base.isequal(a::WignerContribution{C}, b::WignerContribution{C}) where {C<:Number}
    return isequal(a.coefficient, b.coefficient) && isequal(a.kinematic, b.kinematic)
end
Base.:(==)(a::WignerContribution, b::WignerContribution) = isequal(a, b)
function Base.hash(contribution::WignerContribution, h::UInt)
    return hash(
        WignerContribution, hash(contribution.kinematic, hash(contribution.coefficient, h))
    )
end

"""Concrete collection of Wigner-space diagrams at a fixed gradient order/context."""
struct WignerDiagrams{C<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext}
    diagrams::Dict{WignerDiagram{S,E1,E2,G,Nothing,Ctx},Vector{WignerContribution{C}}}
    context::Ctx
end

function WignerDiagrams{C,S,E1,E2,G,Ctx}(
    context::Ctx
) where {C<:Number,S<:Statistics,E1,E2,G,Ctx<:AbstractWignerContext}
    K = WignerDiagram{S,E1,E2,G,Nothing,Ctx}
    return WignerDiagrams{C,S,E1,E2,G,Ctx}(Dict{K,Vector{WignerContribution{C}}}(), context)
end

function SmallCollections.default(
    ::Type{WignerDiagrams{C,S,E1,E2,0,HomogeneousWignerContext}}
) where {C<:Number,S<:Statistics,E1,E2}
    return WignerDiagrams{C,S,E1,E2,0,HomogeneousWignerContext}(HomogeneousWignerContext())
end

Base.length(collection::WignerDiagrams) = length(collection.diagrams)
Base.isempty(collection::WignerDiagrams) = isempty(collection.diagrams)
Base.iszero(collection::WignerDiagrams) = isempty(collection.diagrams)
Base.iterate(collection::WignerDiagrams) = iterate(collection.diagrams)
Base.iterate(collection::WignerDiagrams, state) = iterate(collection.diagrams, state)
wigner_context(collection::WignerDiagrams) = collection.context
function Base.eltype(::Type{WignerDiagrams{C,S,E1,E2,G,Ctx}}) where {C,S,E1,E2,G,Ctx}
    return Pair{WignerDiagram{S,E1,E2,G,Nothing,Ctx},Vector{WignerContribution{C}}}
end
function Base.isequal(
    a::WignerDiagrams{C,S,E1,E2,G,Ctx}, b::WignerDiagrams{C,S,E1,E2,G,Ctx}
) where {C,S,E1,E2,G,Ctx}
    return isequal(a.diagrams, b.diagrams) && isequal(a.context, b.context)
end
Base.:(==)(a::WignerDiagrams, b::WignerDiagrams) = isequal(a, b)
function Base.hash(collection::WignerDiagrams, h::UInt)
    return hash(collection.context, hash(collection.diagrams, h))
end

gradient_order(::WignerDiagrams{C,S,E1,E2,G}) where {C,S,E1,E2,G} = Val(G)

@inline function _homogeneous_wigner_diagram(
    diagram::FourierDiagram{S,E1,E2,K}
)::WignerDiagram{S,E1,E2,0,K,HomogeneousWignerContext} where {S<:Statistics,E1,E2,K}
    external_momentum_count(diagram) == 1 || throw(
        ArgumentError(
            "homogeneous Wigner representation requires exactly one external momentum"
        ),
    )
    basis = momentum_basis(diagram)
    isempty(basis.variables) && error("Wigner input has no external momentum basis")
    return WignerDiagram{S,E1,E2,0,K,HomogeneousWignerContext}(
        diagram, basis[1], HomogeneousWignerContext()
    )
end

function _homogeneous_wigner_diagrams(
    diagrams::FourierDiagrams{C,S,E1,E2}
)::WignerDiagrams{
    C,S,E1,E2,0,HomogeneousWignerContext
} where {C<:Number,S<:Statistics,E1,E2}
    context = HomogeneousWignerContext()
    out = WignerDiagrams{C,S,E1,E2,0,HomogeneousWignerContext}(context)
    for (graph, contributions) in diagrams
        wigner_graph = _homogeneous_wigner_diagram(graph)
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

@noinline function _unsupported_wigner_gradient(::Val{G}) where {G}
    return throw(
        ArgumentError(
            "unsupported Wigner gradient order $G; the initial implementation supports only Val(0)",
        ),
    )
end

_wigner_transform(diagram::FourierDiagram, ::Val{0}) = _homogeneous_wigner_diagram(diagram)
function _wigner_transform(::FourierDiagram, gradient_order::Val{G}) where {G}
    return _unsupported_wigner_gradient(gradient_order)
end

function _wigner_transform(diagrams::FourierDiagrams, ::Val{0})
    return _homogeneous_wigner_diagrams(diagrams)
end
function _wigner_transform(::FourierDiagrams, gradient_order::Val{G}) where {G}
    return _unsupported_wigner_gradient(gradient_order)
end

function wigner_transform(diagram::FourierDiagram; gradient_order::Val=Val(0))
    return _wigner_transform(diagram, gradient_order)
end

function wigner_transform(diagrams::FourierDiagrams; gradient_order::Val=Val(0))
    return _wigner_transform(diagrams, gradient_order)
end

"""Homogeneous Wigner-space dressed propagator at statically visible gradient order `G`."""
struct WignerDressedPropagator{C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    keldysh::WignerDiagrams{C,S,E1,E2,G,Ctx}
    retarded::WignerDiagrams{C,S,E1,E2,G,Ctx}
    advanced::WignerDiagrams{C,S,E1,E2,G,Ctx}
    parameter::ParameterMonomial
    target::FieldFamily{S}
    context::Ctx
end

order(::WignerDressedPropagator{C,S,O}) where {C,S,O} = O
statistics(::WignerDressedPropagator{C,S}) where {C,S} = S
parameters(G::WignerDressedPropagator) = G.parameter
target_family(G::WignerDressedPropagator) = G.target
gradient_order(::WignerDressedPropagator{C,S,O,E1,E2,G}) where {C,S,O,E1,E2,G} = Val(G)
wigner_context(G::WignerDressedPropagator) = G.context
function Base.isequal(
    a::WignerDressedPropagator{C,S,O,E1,E2,G,Ctx},
    b::WignerDressedPropagator{C,S,O,E1,E2,G,Ctx},
) where {C,S,O,E1,E2,G,Ctx}
    return isequal(a.keldysh, b.keldysh) &&
           isequal(a.retarded, b.retarded) &&
           isequal(a.advanced, b.advanced) &&
           isequal(a.parameter, b.parameter) &&
           isequal(a.target, b.target) &&
           isequal(a.context, b.context)
end
Base.:(==)(a::WignerDressedPropagator, b::WignerDressedPropagator) = isequal(a, b)
function Base.hash(G::WignerDressedPropagator, h::UInt)
    return hash((G.keldysh, G.retarded, G.advanced, G.parameter, G.target, G.context), h)
end

function _wigner_transform(
    G::FourierDressedPropagator{C,S,O,E1,E2}, ::Val{0}
)::WignerDressedPropagator{
    C,S,O,E1,E2,0,HomogeneousWignerContext
} where {C<:Number,S<:Statistics,O,E1,E2}
    context = HomogeneousWignerContext()
    return WignerDressedPropagator{C,S,O,E1,E2,0,HomogeneousWignerContext}(
        _homogeneous_wigner_diagrams(G.keldysh),
        _homogeneous_wigner_diagrams(G.retarded),
        _homogeneous_wigner_diagrams(G.advanced),
        G.parameter,
        G.target,
        context,
    )
end
function _wigner_transform(::FourierDressedPropagator, gradient_order::Val{G}) where {G}
    return _unsupported_wigner_gradient(gradient_order)
end

function wigner_transform(G::FourierDressedPropagator; gradient_order::Val=Val(0))
    return _wigner_transform(G, gradient_order)
end

function wigner_transform(G::DressedPropagator; gradient_order::Val=Val(0))
    return wigner_transform(fourier_transform(G); gradient_order=gradient_order)
end

"""Homogeneous Wigner-space self-energy at statically visible gradient order `G`."""
struct WignerSelfEnergy{C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    keldysh::WignerDiagrams{C,S,E1,E2,G,Ctx}
    retarded::WignerDiagrams{C,S,E1,E2,G,Ctx}
    advanced::WignerDiagrams{C,S,E1,E2,G,Ctx}
    parameter::ParameterMonomial
    target::FieldFamily{S}
    context::Ctx
end

order(::WignerSelfEnergy{C,S,O}) where {C,S,O} = O
statistics(::WignerSelfEnergy{C,S}) where {C,S} = S
parameters(Σ::WignerSelfEnergy) = Σ.parameter
target_family(Σ::WignerSelfEnergy) = Σ.target
gradient_order(::WignerSelfEnergy{C,S,O,E1,E2,G}) where {C,S,O,E1,E2,G} = Val(G)
wigner_context(Σ::WignerSelfEnergy) = Σ.context
function Base.isequal(
    a::WignerSelfEnergy{C,S,O,E1,E2,G,Ctx}, b::WignerSelfEnergy{C,S,O,E1,E2,G,Ctx}
) where {C,S,O,E1,E2,G,Ctx}
    return isequal(a.keldysh, b.keldysh) &&
           isequal(a.retarded, b.retarded) &&
           isequal(a.advanced, b.advanced) &&
           isequal(a.parameter, b.parameter) &&
           isequal(a.target, b.target) &&
           isequal(a.context, b.context)
end
Base.:(==)(a::WignerSelfEnergy, b::WignerSelfEnergy) = isequal(a, b)
function Base.hash(Σ::WignerSelfEnergy, h::UInt)
    return hash((Σ.keldysh, Σ.retarded, Σ.advanced, Σ.parameter, Σ.target, Σ.context), h)
end

function _wigner_transform(
    Σ::FourierSelfEnergy{C,S,O,E1,E2}, ::Val{0}
)::WignerSelfEnergy{
    C,S,O,E1,E2,0,HomogeneousWignerContext
} where {C<:Number,S<:Statistics,O,E1,E2}
    context = HomogeneousWignerContext()
    return WignerSelfEnergy{C,S,O,E1,E2,0,HomogeneousWignerContext}(
        _homogeneous_wigner_diagrams(Σ.keldysh),
        _homogeneous_wigner_diagrams(Σ.retarded),
        _homogeneous_wigner_diagrams(Σ.advanced),
        Σ.parameter,
        Σ.target,
        context,
    )
end
function _wigner_transform(::FourierSelfEnergy, gradient_order::Val{G}) where {G}
    return _unsupported_wigner_gradient(gradient_order)
end

function wigner_transform(Σ::FourierSelfEnergy; gradient_order::Val=Val(0))
    return _wigner_transform(Σ, gradient_order)
end

function wigner_transform(::SelfEnergy; gradient_order::Val=Val(0))
    return throw(
        ArgumentError(
            "cannot safely Wigner-transform an amputated coordinate SelfEnergy; transform the DressedPropagator first, construct SelfEnergy(fourier_transform(G)), then call wigner_transform",
        ),
    )
end

function _empty_wigner_diagrams(
    ::Type{WignerDiagrams{C,S,E1,E2,0,HomogeneousWignerContext}}
) where {C<:Number,S<:Statistics,E1,E2}
    return WignerDiagrams{C,S,E1,E2,0,HomogeneousWignerContext}(HomogeneousWignerContext())
end

function matrix(
    G::WignerDressedPropagator{C,Boson,O,E1,E2,0,HomogeneousWignerContext}
) where {C<:Number,O,E1,E2}
    D = WignerDiagrams{C,Boson,E1,E2,0,HomogeneousWignerContext}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = G.keldysh
    result[1, 2] = G.retarded
    result[2, 1] = G.advanced
    result[2, 2] = _empty_wigner_diagrams(D)
    return result
end

function matrix(
    G::WignerDressedPropagator{C,Fermion,O,E1,E2,0,HomogeneousWignerContext}
) where {C<:Number,O,E1,E2}
    D = WignerDiagrams{C,Fermion,E1,E2,0,HomogeneousWignerContext}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = G.retarded
    result[1, 2] = G.keldysh
    result[2, 1] = _empty_wigner_diagrams(D)
    result[2, 2] = G.advanced
    return result
end

function matrix(
    Σ::WignerSelfEnergy{C,Boson,O,E1,E2,0,HomogeneousWignerContext}
) where {C<:Number,O,E1,E2}
    D = WignerDiagrams{C,Boson,E1,E2,0,HomogeneousWignerContext}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = _empty_wigner_diagrams(D)
    result[1, 2] = Σ.advanced
    result[2, 1] = Σ.retarded
    result[2, 2] = Σ.keldysh
    return result
end

function matrix(
    Σ::WignerSelfEnergy{C,Fermion,O,E1,E2,0,HomogeneousWignerContext}
) where {C<:Number,O,E1,E2}
    D = WignerDiagrams{C,Fermion,E1,E2,0,HomogeneousWignerContext}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = Σ.retarded
    result[1, 2] = Σ.keldysh
    result[2, 1] = _empty_wigner_diagrams(D)
    result[2, 2] = Σ.advanced
    return result
end
