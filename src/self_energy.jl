"Compute the self-energy component from the external propagator types."
function self_energy_type(::Type{Boson}, dict::SmallCollections.SmallDict)
    if is_keldysh(dict[:out]) && is_advanced(dict[:in])
        return PropagatorType.Advanced
    elseif is_retarded(dict[:out]) && is_keldysh(dict[:in])
        return PropagatorType.Retarded
    elseif is_retarded(dict[:out]) && is_advanced(dict[:in])
        return PropagatorType.Keldysh
    else
        error("Classical-Classical for self-energy should be zero.")
    end
end

"Construct the self-energy from irreducible diagrams."
function construct_self_energy!(
    self_energy::SmallCollections.SmallDict,
    diagrams::Diagrams{C,S,E,E2},
) where {C<:Number,S<:Statistics,E,E2}
    for (diagram, prefactor) in diagrams
        _contractions = contractions(diagram)

        if !is_irreducible(_contractions)
            continue
        end

        positions = position_category.(_contractions)
        types_p = propagator_type.(_contractions)
        dict = SmallCollections.SmallDict{E,Symbol,PropagatorType.T}(
            p => t for (p, t) in zip(positions, types_p)
        )

        if is_keldysh(dict[:out]) && is_keldysh(dict[:in])
            continue
        end

        bulk_propagators = FixedVector{E - 2,Edge{S}}(
            edge for edge in _contractions if is_bulk(edge)
        )
        push!(
            self_energy[self_energy_type(S, dict)],
            Diagram(bulk_propagators, Val(E2)),
            prefactor,
        )
    end
    return self_energy
end

"""
$(DocStringExtensions.TYPEDEF)

Self-energy components in the Retarded-Advanced-Keldysh basis. Coefficient representation,
statistics, perturbation order, and diagram shape are encoded in the type.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SelfEnergy{C<:Number,S<:Statistics,O,E1,E2}
    "The Keldysh component of the self-energy."
    keldysh::Diagrams{C,S,E1,E2}
    "The retarded component of the self-energy."
    retarded::Diagrams{C,S,E1,E2}
    "The advanced component of the self-energy."
    advanced::Diagrams{C,S,E1,E2}
    "Parameters of the perturbation series"
    parameter::CSym
end

order(::SelfEnergy{C,S,O}) where {C,S,O} = O
statistics(::SelfEnergy{C,S}) where {C,S} = S
parameters(Σ::SelfEnergy) = Σ.parameter

function self_energy_result_type(
    ::Type{DressedPropagator{C,S,O,E1,E2}}
) where {C<:Number,S<:Statistics,O,E1,E2}
    return SelfEnergy{C,S,O,E1 - 2,max_edges(O)}
end

function _self_energy(
    G::DressedPropagator{C,Boson,O,E1,E2}
) where {C<:Number,O,E1,E2}
    SE = E1 - 2
    ST = max_edges(O)
    D = Diagrams{C,Boson,SE,ST}
    self_energy = SmallCollections.SmallDict{3,PropagatorType.T,D}((
        PropagatorType.Advanced => D(),
        PropagatorType.Retarded => D(),
        PropagatorType.Keldysh => D(),
    ))
    construct_self_energy!(self_energy, G.keldysh)

    _simplify_prefactors!(self_energy[PropagatorType.Keldysh])
    _simplify_prefactors!(self_energy[PropagatorType.Retarded])
    _simplify_prefactors!(self_energy[PropagatorType.Advanced])
    return SelfEnergy{C,Boson,O,SE,ST}(
        self_energy[PropagatorType.Keldysh],
        self_energy[PropagatorType.Retarded],
        self_energy[PropagatorType.Advanced],
        G.parameter,
    )
end

SelfEnergy(G::DressedPropagator) = _self_energy(G)

# Transitional call shape while repository call sites move to `SelfEnergy(G)`.
function SelfEnergy(G::DressedPropagator{C,S,O}, ::Val{O}) where {C,S,O}
    return SelfEnergy(G)
end

"""
    matrix(Σ::SelfEnergy)

Return the bosonic Retarded-Advanced-Keldysh self-energy matrix
```math
\\hat{\\Sigma}=\\begin{pmatrix}0&\\Sigma^A\\\\\\Sigma^R&\\Sigma^K\\end{pmatrix}.
```
"""
function matrix(Σ::SelfEnergy{C,Boson,O,E1,E2}) where {C<:Number,O,E1,E2}
    return matrix(Boson, Σ)
end
function matrix(
    ::Type{Boson}, Σ::SelfEnergy{C,Boson,O,E1,E2}
) where {C<:Number,O,E1,E2}
    D = Diagrams{C,Boson,E1,E2}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = structural_zero(Boson, D)
    result[1, 2] = Σ.advanced
    result[2, 1] = Σ.retarded
    result[2, 2] = Σ.keldysh
    return result
end

"""Collection of self-energies with distinct perturbation-parameter monomials."""
struct SelfEnergySum{K,ΣT,O}
    arguments::Dict{K,ΣT}
end

SymbolicUtils.arguments(d::SelfEnergySum) = d.arguments
order(::SelfEnergySum{K,ΣT,O}) where {K,ΣT,O} = O
parameters(d::SelfEnergySum) = map(Σ -> Σ.parameter, values(arguments(d)))

function SelfEnergy(G::DressedPropagatorSum{K,GS,O}) where {K,GS,O}
    ΣT = self_energy_result_type(GS)
    dict = Dict{K,ΣT}(key => SelfEnergy(val) for (key, val) in arguments(G))
    return SelfEnergySum{K,ΣT,O}(dict)
end

# Transitional call shape while repository call sites move to `SelfEnergy(G)`.
function SelfEnergy(G::DressedPropagatorSum{K,GS,O}, ::Val{O}) where {K,GS,O}
    return SelfEnergy(G)
end
