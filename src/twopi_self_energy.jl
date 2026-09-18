"""Formal two-point diagram obtained by cutting one line of a 2PI vacuum skeleton."""
struct TwoPISelfEnergyDiagram{S<:Statistics,E,E2}
    contractions::FixedVector{E,Contraction{S}}
    row::KeldyshIndex.T
    column::KeldyshIndex.T
    topology::FixedVector{E2,Int}
end

function Base.isequal(
    a::TwoPISelfEnergyDiagram{S,E,E2}, b::TwoPISelfEnergyDiagram{S,E,E2}
) where {S,E,E2}
    return a.row === b.row &&
           a.column === b.column &&
           isequal(a.contractions, b.contractions)
end
function Base.:(==)(
    a::TwoPISelfEnergyDiagram{S,E,E2}, b::TwoPISelfEnergyDiagram{S,E,E2}
) where {S,E,E2}
    return isequal(a, b)
end
function Base.hash(d::TwoPISelfEnergyDiagram, h::UInt)
    return hash(
        TwoPISelfEnergyDiagram, hash(d.column, hash(d.row, hash(d.contractions, h)))
    )
end

twopi_self_energy_contractions(d::TwoPISelfEnergyDiagram) = d.contractions
twopi_self_energy_component(d::TwoPISelfEnergyDiagram) = (d.row, d.column)

"""
$(DocStringExtensions.TYPEDEF)

Formal proper self-energy obtained by graph-functional differentiation of a
[`TwoPIEffectiveAction`](@ref).

The representation still lives on the formal full-propagator space: internal quantum--quantum
contractions are retained until an explicit physical projection is requested. The two Keldysh
indices stored by each diagram are the row and column indices of the self-energy after the
transpose in ``Σ_{a b} ∝ δΓ₂/δG_{b a}``.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoPISelfEnergyFunctional{C<:Number,S<:Statistics,O,E,E2}
    "Canonical open two-point diagrams with exact derivative coefficients."
    diagrams::Dict{TwoPISelfEnergyDiagram{S,E,E2},C}
    "Canonical perturbation-parameter monomial."
    parameter::ParameterMonomial
    "Field family with respect to whose full propagator Γ₂ was differentiated."
    target::FieldFamily{S}
end

order(::TwoPISelfEnergyFunctional{C,S,O}) where {C,S,O} = O
statistics(::TwoPISelfEnergyFunctional{C,S}) where {C,S} = S
parameters(Σ::TwoPISelfEnergyFunctional) = Σ.parameter
target_family(Σ::TwoPISelfEnergyFunctional) = Σ.target
Base.length(Σ::TwoPISelfEnergyFunctional) = length(Σ.diagrams)
Base.iszero(Σ::TwoPISelfEnergyFunctional) = isempty(Σ.diagrams)
Base.iterate(Σ::TwoPISelfEnergyFunctional) = iterate(Σ.diagrams)
Base.iterate(Σ::TwoPISelfEnergyFunctional, state) = iterate(Σ.diagrams, state)

"""Return a copy of the exact formal 2PI self-energy derivative terms."""
twopi_self_energy_terms(Σ::TwoPISelfEnergyFunctional) = copy(Σ.diagrams)

function _twopi_cut_anchors(cut::Contraction{S}) where {S<:Statistics}
    out_vertex = position(cut.out)
    in_vertex = position(cut.in)
    out_anchor = Contraction(cut.out(Out()), bar(cut.out)(out_vertex))
    in_anchor = Contraction(bar(cut.in)(in_vertex), cut.in(In()))
    return out_anchor, in_anchor
end

function _twopi_canonical_cut(
    contractions::FixedVector{E,Contraction{S}}, cut_index::Int, ::Val{E2}
) where {S<:Statistics,E,E2}
    cut = contractions[cut_index]
    out_anchor, in_anchor = _twopi_cut_anchors(cut)

    anchored = Contraction{S}[]
    sizehint!(anchored, E + 1)
    for index in eachindex(contractions)
        index == cut_index && continue
        push!(anchored, contractions[index])
    end
    push!(anchored, out_anchor)
    push!(anchored, in_anchor)

    canonical = canonicalize(anchored)
    internal = Contraction{S}[
        contraction for contraction in canonical if is_bulk(contraction)
    ]
    length(internal) == E - 1 || error("2PI line cut produced an invalid anchored graph")
    key = sorted_wick_key(internal, Val(E - 1))
    topology = legacy_topology(collect(key), Val(E2))
    row = keldysh_index(cut.in)
    column = keldysh_index(cut.out)
    return TwoPISelfEnergyDiagram{S,E - 1,E2}(key, row, column, topology)
end

# KC represents a complex bosonic propagator as one oriented ψ -> ψ̄ line. In this
# non-Nambu representation the 2PI stationary functional contains i Tr rather than
# (i/2) Tr over a doubled field space, hence Σ = i δΓ₂/δGᵀ. A Nambu-doubled
# representation instead carries the familiar 2i prefactor.
_twopi_derivative_factor(::Type{Boson}, ::Type{C}) where {C<:Number} = convert(C, im)
function _twopi_derivative_factor(::Type{Fermion}, ::Type{C}) where {C<:Number}
    return throw(
        ArgumentError(
            "fermionic 2PI functional-derivative normalization is not yet certified for KC's oriented propagator representation",
        ),
    )
end

"""
    twopi_self_energy(Γ₂, target)

Differentiate a 2PI effective action with respect to the formal full propagator of `target`.
For KC's oriented complex-boson propagator representation this implements

```math
Σ_{a b}(x,x') = i δΓ₂/δG_{b a}(x',x).
```

This is the non-Nambu form of the 2PI functional derivative. If a complex field is instead written
as an explicitly doubled Nambu matrix, the effective-action trace carries an additional factor
`1/2` and the equivalent convention is `Σ = 2i δΓ₂/δGᵀ`.

Every target-family propagator line is cut once. The two cut vertices are temporarily attached to
`Out()` and `In()` marker legs before canonicalization, so open two-point graphs are compared with
both external attachment points distinguished. The markers are stripped only after canonical
relabeling. Equivalent cuts therefore combine with their exact line multiplicity.

This operation does not impose `G_qq = 0` on the remaining internal lines. Use
[`physical_self_energy`](@ref) for the subsequent physical Keldysh projection.
"""
function twopi_self_energy(
    Γ::TwoPIEffectiveAction{C,S,O,E,E2}, target::FieldFamily{S}
) where {C<:Number,S<:Statistics,O,E,E2}
    target in field_families(Γ) ||
        throw(ArgumentError("target field family is absent from the 2PI effective action"))

    derivative_factor = _twopi_derivative_factor(S, C)
    diagrams = Dict{TwoPISelfEnergyDiagram{S,E - 1,E2},C}()
    for (vacuum, coefficient) in Γ
        contractions = twopi_contractions(vacuum)
        for cut_index in eachindex(contractions)
            cut = contractions[cut_index]
            isequal(field_family(cut.out), target) || continue
            open_diagram = _twopi_canonical_cut(contractions, cut_index, Val(E2))
            contribution = _simplify(derivative_factor * coefficient)
            combined = _simplify(get(diagrams, open_diagram, zero(C)) + contribution)
            if iszero(combined)
                delete!(diagrams, open_diagram)
            else
                diagrams[open_diagram] = combined
            end
        end
    end
    return TwoPISelfEnergyFunctional{C,S,O,E - 1,E2}(diagrams, parameters(Γ), target)
end

function _twopi_physical_component(row::KeldyshIndex.T, column::KeldyshIndex.T)::UInt8
    if row === Classical && column === Classical
        return 0x00
    elseif row === Quantum && column === Classical
        return 0x01
    elseif row === Classical && column === Quantum
        return 0x02
    elseif row === Quantum && column === Quantum
        return 0x03
    end
    return error("unknown Keldysh self-energy component")
end

function _twopi_internal_is_physical(contractions)
    all(contraction_filter, contractions) || return false
    all(regular, contractions) || return false
    return !has_zero_loop(collect(contractions))
end

"""
    physical_self_energy(Σ₂PI)

Restrict a formal bosonic 2PI self-energy derivative to the physical Keldysh propagator manifold.
The structural ``Σ_cl,cl`` component is discarded, residual internal quantum--quantum propagators
and causal zero loops vanish, and surviving components are converted to the existing
[`SelfEnergy`](@ref) representation.

Physical projection is deliberately separate from [`twopi_self_energy`](@ref): cutting a formal
`G_qq` line is what generates the physical Keldysh self-energy ``Σ_qq``.
"""
function physical_self_energy(
    Σ::TwoPISelfEnergyFunctional{C,Boson,O,E,E2}
) where {C<:Number,O,E,E2}
    D = Diagrams{C,Boson,E,E2}
    components = SmallCollections.SmallDict{3,PropagatorType.T,D}((
        PropagatorType.Advanced => D(),
        PropagatorType.Retarded => D(),
        PropagatorType.Keldysh => D(),
    ))

    for (formal_diagram, coefficient) in Σ
        component = _twopi_physical_component(formal_diagram.row, formal_diagram.column)
        iszero(component) && continue
        internal = formal_diagram.contractions
        _twopi_internal_is_physical(internal) || continue
        E > 0 || throw(ArgumentError("cannot project a zero-line 2PI self-energy diagram"))
        diagram = Diagram(collect(internal), Val(E), Val(E2))
        if component == 0x01
            push!(components[PropagatorType.Retarded], diagram, coefficient)
        elseif component == 0x02
            push!(components[PropagatorType.Advanced], diagram, coefficient)
        else
            push!(components[PropagatorType.Keldysh], diagram, coefficient)
        end
    end

    _simplify_prefactors!(components[PropagatorType.Keldysh])
    _simplify_prefactors!(components[PropagatorType.Retarded])
    _simplify_prefactors!(components[PropagatorType.Advanced])
    return SelfEnergy{C,Boson,O,E,E2}(
        components[PropagatorType.Keldysh],
        components[PropagatorType.Retarded],
        components[PropagatorType.Advanced],
        parameters(Σ),
        target_family(Σ),
    )
end

"""
    SelfEnergy(Γ₂::TwoPIEffectiveAction, target::FieldFamily)

Lower a bosonic 2PI effective action to KC's ordinary physical [`SelfEnergy`](@ref)
representation for `target`.

This is the high-level 2PI counterpart of `SelfEnergy(G::DressedPropagator)`. It is deliberately
a thin semantic lowering through [`twopi_self_energy`](@ref) followed by
[`physical_self_energy`](@ref); no separate graph or projection path is introduced.
"""
function SelfEnergy(
    Γ::TwoPIEffectiveAction{C,Boson,O,E,E2}, target::FieldFamily{Boson}
) where {C<:Number,O,E,E2}
    return physical_self_energy(twopi_self_energy(Γ, target))
end

"""
    SelfEnergy(Γ₂::TwoPIEffectiveAction)

Lower a single-family bosonic 2PI effective action to the ordinary physical
[`SelfEnergy`](@ref) representation.

For a multi-family effective action, such as the Hubbard--Stratonovich `G²D` theory, the target
family is physically meaningful and must be supplied explicitly with `SelfEnergy(Γ₂, target)`.
"""
function SelfEnergy(Γ::TwoPIEffectiveAction{C,Boson,O,E,E2}) where {C<:Number,O,E,E2}
    families = field_families(Γ)
    length(families) == 1 || throw(
        ArgumentError(
            "a multi-family 2PI effective action requires an explicit self-energy target",
        ),
    )
    return SelfEnergy(Γ, only(families))
end
