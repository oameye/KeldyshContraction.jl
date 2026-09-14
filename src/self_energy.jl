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
    self_energy::SmallCollections.SmallDict, diagrams::Diagrams{C,S,E,E2}
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

Fixed-order one-particle-irreducible two-point self-energy in the
Retarded-Advanced-Keldysh basis. The result stores semantic Keldysh, retarded, and advanced
components and preserves the perturbative parameter and external physical field family of the
input propagator.

Use [`keldysh_component`](@ref), [`retarded_component`](@ref), and
[`advanced_component`](@ref) for component access. `SelfEnergy` constructs the explicit
perturbative 1PI kernel; it does not solve a self-consistent Dyson equation.

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
    "Canonical perturbation-parameter monomial"
    parameter::ParameterMonomial
    "Physical field family of the amputated external two-point function"
    target::FieldFamily{S}
end

order(::SelfEnergy{C,S,O}) where {C,S,O} = O
statistics(::SelfEnergy{C,S}) where {C,S} = S
parameters(Σ::SelfEnergy) = Σ.parameter
target_family(Σ::SelfEnergy) = Σ.target

function self_energy_result_type(
    ::Type{DressedPropagator{C,S,O,E1,E2}}
) where {C<:Number,S<:Statistics,O,E1,E2}
    return SelfEnergy{C,S,O,E1 - 2,max_edges(O)}
end

function _self_energy(G::DressedPropagator{C,Boson,O,E1,E2}) where {C<:Number,O,E1,E2}
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
        G.target,
    )
end

"""
    SelfEnergy(G::DressedPropagator)

Extract the fixed-order 1PI self-energy from the perturbative dressed two-point function `G`.
The two external propagator lines are amputated, reducible two-point diagrams are discarded,
and the exact diagrammatic prefactors and statistics are retained.
"""
SelfEnergy(G::DressedPropagator) = _self_energy(G)

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
function matrix(::Type{Boson}, Σ::SelfEnergy{C,Boson,O,E1,E2}) where {C<:Number,O,E1,E2}
    D = Diagrams{C,Boson,E1,E2}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = structural_zero(Boson, D)
    result[1, 2] = Σ.advanced
    result[2, 1] = Σ.retarded
    result[2, 2] = Σ.keldysh
    return result
end

"""Collection of self-energies with distinct perturbation-parameter monomials."""
struct SelfEnergySum{ΣT,O}
    arguments::Dict{ParameterMonomial,ΣT}
end

SymbolicUtils.arguments(d::SelfEnergySum) = d.arguments
Base.getindex(d::SelfEnergySum, parameter) = d.arguments[parameter_monomial(parameter)]
order(::SelfEnergySum{ΣT,O}) where {ΣT,O} = O
parameters(d::SelfEnergySum) = collect(keys(d.arguments))

"""
    SelfEnergy(G::DressedPropagatorSum)

Extract the 1PI self-energy independently for every perturbative parameter sector in `G` and
preserve the sector keys for downstream Fourier, Wigner, and kinetic transformations.
"""
function SelfEnergy(G::DressedPropagatorSum{GS,O}) where {GS,O}
    ΣT = self_energy_result_type(GS)
    dict = Dict{ParameterMonomial,ΣT}(key => SelfEnergy(val) for (key, val) in arguments(G))
    return SelfEnergySum{ΣT,O}(dict)
end
