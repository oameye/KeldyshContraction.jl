"""
Return the canonical endpoint regularisations for one relative equal-time shift.

The canonical representative preserves the exact difference `ε_out-ε_in` while choosing a
unique spelling. Relative shifts outside `-2:2` cannot arise from the package's three-valued
endpoint regularisation and are rejected explicitly.
"""
@inline function _canonical_equal_time_regularisations(relative::Int)
    if relative == 2
        return (Regularisation.Plus, Regularisation.Minus)
    elseif relative == 1
        return (Regularisation.Plus, Regularisation.Zero)
    elseif iszero(relative)
        return (Regularisation.Zero, Regularisation.Zero)
    elseif relative == -1
        return (Regularisation.Minus, Regularisation.Zero)
    elseif relative == -2
        return (Regularisation.Minus, Regularisation.Plus)
    end
    throw(
        ArgumentError("equal-time relative regularisation must lie in -2:2, got $relative")
    )
end

@inline function _uses_equal_time_causal_regularisation(edge::Edge)
    same_position = isequal(position(edge.out), position(edge.in))
    causal = is_retarded(edge) || is_advanced(edge)
    return same_position && causal
end

"""
Canonicalize only the spelling of an equal-position retarded/advanced edge.

This transform is deliberately separate from `Edge.isequal`/`hash`: diagram generation must
retain endpoint Trotter provenance and must not merge Wick contributions merely because two
causal edges have the same relative shift. The canonical form is intended for physical
comparison after diagram coefficients have already been assembled.
"""
function canonicalize_equal_time_regularisation(edge::Edge{S}) where {S<:Statistics}
    _uses_equal_time_causal_regularisation(edge) || return edge
    relative = subtraction(regularisations(edge))
    out_reg, in_reg = _canonical_equal_time_regularisations(relative)
    out = reconstruct(edge.out; regularisation=out_reg)
    incoming = reconstruct(edge.in; regularisation=in_reg)
    return Edge{S}(out, incoming, edge.edgetype, edge.momenta)
end

function canonicalize_equal_time_regularisation(
    diagram::Diagram{S,E1,E2}
) where {S<:Statistics,E1,E2}
    edges = FixedVector{E1,Edge{S}}(
        canonicalize_equal_time_regularisation(edge) for edge in contractions(diagram)
    )
    return Diagram{S,E1,E2}(edges, topology(diagram))
end

function _canonicalize_equal_time_causal_diagram(
    diagram::Diagram{S,E1,E2}
) where {S<:Statistics,E1,E2}
    edges = Vector{Edge{S}}(undef, E1)
    sign = 1
    for (index, edge) in enumerate(contractions(diagram))
        if _uses_equal_time_causal_regularisation(edge) && is_advanced(edge)
            edges[index] = canonicalize_equal_time_regularisation(adjoint(edge))
            sign = -sign
        else
            edges[index] = canonicalize_equal_time_regularisation(edge)
        end
    end
    fixed_edges = FixedVector{E1,Edge{S}}(edge for edge in edges)
    return Diagram{S,E1,E2}(fixed_edges, topology(diagram)), sign
end

"""
Canonicalize equal-time causal structure in a completed diagram collection.

Endpoint spellings are first reduced to their unique relative-shift representatives. An
equal-position advanced line is then represented by the retarded line using
`Gᴬ(y,y) = -Gᴿ(y,y)`; the minus sign is carried by the completed diagram coefficient rather
than hidden in `Edge` identity. Equivalent keys are combined only in the returned copy, after
Wick/topology multiplicities have already been assembled. The input collection is unchanged.
"""
function canonicalize_equal_time_regularisation(
    diagrams::Diagrams{C,S,E1,E2}
) where {C<:Number,S<:Statistics,E1,E2}
    out = Diagrams{C,S,E1,E2}()
    for (diagram, coefficient) in diagrams
        canonical, sign = _canonicalize_equal_time_causal_diagram(diagram)
        push!(out, canonical, sign * coefficient)
    end
    return out
end
