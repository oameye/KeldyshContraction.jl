function _skeleton_fourier_diagrams(
    diagrams::FourierDiagrams{C,S,E1,E2}
) where {C<:Number,S<:Statistics,E1,E2}
    out = FourierDiagrams{C,S,E1,E2}()
    for (graph, contributions) in diagrams
        coordinate_edges = contractions(graph.coordinate)
        is_two_particle_irreducible(coordinate_edges) || continue
        out.diagrams[graph] = copy(contributions)
    end
    return out
end

"""
    skeleton_self_energy(G::FourierDressedPropagator)

Extract the 2PI/skeleton subset of the ordinary Fourier-space 1PI self-energy. The graph filter is
applied to the amputated coordinate graph carried by each exact Fourier diagram; momentum routing,
derivative-generated kinematic factors, and diagram coefficients are left untouched.
"""
function skeleton_self_energy(G::FourierDressedPropagator)
    Σ = SelfEnergy(G)
    return typeof(Σ)(
        _skeleton_fourier_diagrams(Σ.keldysh),
        _skeleton_fourier_diagrams(Σ.retarded),
        _skeleton_fourier_diagrams(Σ.advanced),
        Σ.parameter,
        Σ.target,
    )
end
