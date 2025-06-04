const PositionPropagatorType = OrderedCollections.LittleDict{
    AbstractPosition,
    PropagatorType,
    Tuple{<:AbstractPosition,<:AbstractPosition,<:AbstractPosition},
    Tuple{PropagatorType,PropagatorType,PropagatorType},
}

"compute the self-energy type from positions save in `dict`."
function self_energy_type(dict::OrderedCollections.LittleDict)
    # G1K = GK[x1,y1] ΣA[y1,y1] GA[y1,x2] +GR[x1,y1] ΣK GA[y1,x2] + GR[x1,y1]ΣR[y1,y1] GK[y1,x2]
    if is_keldysh(dict[Out()]) && is_advanced(dict[In()])
        return Advanced
    elseif is_retarded(dict[Out()]) && is_keldysh(dict[In()])
        return Retarded
    elseif is_retarded(dict[Out()]) && is_advanced(dict[In()])
        return Keldysh
    else
        @show dict
        error("Classical-Classical for self-energy should be zero.")
    end
end

"Construct the self-energy from diagrams and save in LittleDict self_energy."
function construct_self_energy!(
    self_energy::OrderedCollections.LittleDict, diagrams::Diagrams; order::Int=1
)
    if order > 2
        error("Higher than second order in self-energy is not supported.")
    end

    for (diagram, prefactor) in diagrams
        _contractions = contractions(diagram)

        mult = topology(diagram)
        if !isempty(mult) && first(mult) < order
            continue
        end

        positions = position.(_contractions)
        types_p = propagator_type.(_contractions)
        dict = OrderedCollections.freeze(
            OrderedCollections.OrderedDict(zip(positions, types_p))
        ) # TODO Replace with SmallCollections

        # Find all bulk propagators (edges where both fields are bulk)
        bulk_idxs = findall(isbulk, _contractions)
        bulk_propagators = _contractions[bulk_idxs]
        push!(self_energy[self_energy_type(dict)], Diagram(bulk_propagators), prefactor)
    end
    return self_energy
end

"""
$(TYPEDEF)

A struct representing the self-energy components in the Retarded-Advance-Keldysh basis ([`PropagatorType`](@ref)).
The self-energy is divided into three components: Keldysh, retarded, and advanced.

# Fields
$(FIELDS)
where it assumed that the fields are of type `Union{SymbolicUtils.Symbolic{<:Number}, Number}`.


# Constructor
$(TYPEDSIGNATURES)
Constructs a `SelfEnergy` object from a [`DressedPropagator`](@ref).
The self-energy is computed based on the Keldysh Green's function (`G.keldysh`) and expanded into
its quantum-quantum (`qq`), classical-quantum (`cq`), and quantum-classical (`qc`) components.

"""
struct SelfEnergy{E,T}
    "The Keldysh component of the self-energy."
    keldysh::Diagrams{E,T}
    " The retarded component of the self-energy."
    retarded::Diagrams{E,T}
    "The advanced component of the self-energy."
    advanced::Diagrams{E,T}
end
function SelfEnergy(G::DressedPropagator{E}; order=1) where {E}
    self_energy = OrderedCollections.LittleDict{PropagatorType,Diagrams}((
        Advanced => Diagrams{E-2}(), Retarded => Diagrams{E-2}(), Keldysh => Diagrams{E-2}()
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
    _simplify_prefactors!(self_energy[Keldysh])
    _simplify_prefactors!(self_energy[Retarded])
    _simplify_prefactors!(self_energy[Advanced])
    return SelfEnergy(self_energy[Keldysh], self_energy[Retarded], self_energy[Advanced])
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
    Diagrams[Diagrams{E}() Σ.advanced; Σ.retarded Σ.keldysh]
end
