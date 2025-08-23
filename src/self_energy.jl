"compute the self-energy type from positions save in `dict`."
function self_energy_type(dict::SmallCollections.SmallDict)
    # G1K = GK[x1,y1] ΣA[y1,y1] GA[y1,x2] +GR[x1,y1] ΣK GA[y1,x2] + GR[x1,y1]ΣR[y1,y1] GK[y1,x2]
    if is_keldysh(dict[:out]) && is_advanced(dict[:in])
        return PropagatorType.Advanced
    elseif is_retarded(dict[:out]) && is_keldysh(dict[:in])
        return PropagatorType.Retarded
    elseif is_retarded(dict[:out]) && is_advanced(dict[:in])
        return PropagatorType.Keldysh
    else
        @show dict
        error("Classical-Classical for self-energy should be zero.")
    end
end

"Construct the self-energy from diagrams and save in LittleDict self_energy."
function construct_self_energy!(
    self_energy::SmallCollections.SmallDict, diagrams::Diagrams{E}; order::Int=1
) where {E}
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

        # Find all bulk propagators (edges where both fields are bulk)
        bulk_idxs = findall(is_bulk, _contractions)
        bulk_propagators = _contractions[bulk_idxs]
        push!(self_energy[self_energy_type(dict)], Diagram(bulk_propagators), prefactor)
    end
    return self_energy
end

"""
$(DocStringExtensions.TYPEDEF)

A struct representing the self-energy components in the Retarded-Advance-Keldysh basis ([`PropagatorType`](@ref)).
The self-energy is divided into three components: Keldysh, retarded, and advanced.

# Fields
$(DocStringExtensions.FIELDS)
where it assumed that the fields are of type `Union{SymbolicUtils.Symbolic{<:Number}, Number}`.


# Constructor
$(DocStringExtensions.TYPEDSIGNATURES)
Constructs a `SelfEnergy` object from a [`DressedPropagator`](@ref).
The self-energy is computed based on the Keldysh Green's function (`G.keldysh`) and expanded into
its quantum-quantum (`qq`), classical-quantum (`cq`), and quantum-classical (`qc`) components.

"""
struct SelfEnergy{E1,E2}
    "The Keldysh component of the self-energy."
    keldysh::Diagrams{E1,E2}
    " The retarded component of the self-energy."
    retarded::Diagrams{E1,E2}
    "The advanced component of the self-energy."
    advanced::Diagrams{E1,E2}
    "The order of the self-energy in the perturbation series"
    order::Int64
    "Parameters of the perturbation series"
    parameter::SymbolicUtils.BasicSymbolic{Number}
end
function SelfEnergy(G::DressedPropagator{E}, order=G.order) where {E}
    self_energy = SmallCollections.SmallDict{3,PropagatorType.T,Diagrams}((
        PropagatorType.Advanced => Diagrams{E - 2,topology_length(E)}(),
        PropagatorType.Retarded => Diagrams{E - 2,topology_length(E)}(),
        PropagatorType.Keldysh => Diagrams{E - 2,topology_length(E)}(),
    ))
    construct_self_energy!(self_energy, G.keldysh; order)
    # ^ keldysh GF should contain everything
    # construct_self_energy!(self_energy, G.advanced)
    # construct_self_energy!(self_energy, G.retarded)

    # quantum-quantum is the keldysh term in the self-energy
    # classical-classical is zero

    # G_R(1) = G₀_R Σ_R G₀_R
    # G_A(1) = G₀_A Σ_A G₀_A
    # G_K(1) = G₀_K(x1) Σ_A(y) G₀_A(x2) + G_A(x2) Σ_A(y) G_R(x1) + G_R(x1) Σ_R(y) G_K(x2)
    # G₀_K Σ_K G₀_K = 0
    _simplify_prefactors!(self_energy[PropagatorType.Keldysh])
    _simplify_prefactors!(self_energy[PropagatorType.Retarded])
    _simplify_prefactors!(self_energy[PropagatorType.Advanced])
    return SelfEnergy(
        self_energy[PropagatorType.Keldysh],
        self_energy[PropagatorType.Retarded],
        self_energy[PropagatorType.Advanced],
        order,
        G.parameter,
    )
end

"""
    matrix(Σ::SelfEnergy)

Returns the matrix representation of the self energy `Σ`
in the Retarded-Advanced-Keldysh basis.
```math
\\hat{\\Sigma}\\left(y_1, y_2\\right)=
\\left(\\begin{array}{cc}0 & \\Sigma^A\\left(y_1, y_2\\right) \\\\
\\Sigma^R\\left(y_1, y_2\\right) & \\Sigma^K\\left(y_1, y_2\\right)
\\end{array}
\\right)
```
"""
function matrix(Σ::SelfEnergy{E}) where {E}
    return Diagrams[Diagrams{E,E}() Σ.advanced; Σ.retarded Σ.keldysh]
end

parameters(Σ::SelfEnergy) = Σ.parameter

"""
$(DocStringExtensions.TYPEDEF)

A structure representing a collection of dressed propagator each involving a
different parameter of the expansion.

# Fields
$(DocStringExtensions.FIELDS)

"""
struct SelfEnergySum{Σs}
    "The arguments of the dressed propagator sum. Each involving a different parameter."
    arguments::Dict{SymbolicUtils.BasicSymbolic{Number},Σs}
    "The order of the dressed propagator in the perturbation series"
    order::Int64
end

SymbolicUtils.arguments(d::SelfEnergySum) = d.arguments
order(d::SelfEnergySum) = d.order
parameters(d::SelfEnergySum) = map(G -> G.parameter, arguments(d))

function SelfEnergy(G::DressedPropagatorSum, order=order(G))
    dict = Dict(key => SelfEnergy(val, order) for (key, val) in arguments(G))
    return SelfEnergySum(dict, order)
end
