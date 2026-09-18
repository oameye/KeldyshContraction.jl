"""
    LoopMomentumQuotientWorkspace(expression::OccupationReducedExpression)

Reusable scratch storage for repeated exact dummy-loop quotienting.

The workspace is prepared at the graph and loop capacities required by `expression`. Reusing it
keeps GraphCombinations canonicalization scratch outside the hot quotient call while preserving
exactly the same loop-gauge semantics and output representation as `quotient_loop_momenta`.
"""
struct LoopMomentumQuotientWorkspace
    storage::_LoopGCQuotientWorkspace
end

function LoopMomentumQuotientWorkspace(expression::OccupationReducedExpression)
    graph_capacity, loop_capacity = _loop_gc_capacities(expression)
    return LoopMomentumQuotientWorkspace(
        _LoopGCQuotientWorkspace(graph_capacity, loop_capacity)
    )
end

function _quotient_loop_momenta_graphcombinations(
    expression::OccupationReducedExpression{C,S,O,G,Ctx},
    workspace::_LoopGCQuotientWorkspace,
) where {C<:Number,S<:Statistics,O,G,Ctx<:AbstractWignerContext}
    D = promote_type(C, ComplexRationals, Rational{Int})
    out = Dict{CollisionKernelSector{S},OccupationPolynomial{D,S}}()

    for (sector, polynomial) in occupation_reduced_terms(expression)
        for (occupation_monomial, occupation_coefficient) in polynomial
            for (kinematic_monomial, kinematic_coefficient) in kinematic_factor(sector)
                atom_sector = _kinematic_atom_sector(sector, kinematic_monomial)
                external_index = _prepare_matrixfree_gc_gauge!(
                    atom_sector, occupation_monomial, workspace
                )
                transformed_sector, support_factor = _matrixfree_gc_sector(
                    atom_sector, workspace, external_index
                )
                transformed_monomial = _matrixfree_gc_occupation_monomial(
                    occupation_monomial, workspace, external_index
                )
                transformed_coefficient =
                    convert(D, occupation_coefficient) *
                    convert(D, kinematic_coefficient) *
                    convert(D, support_factor)
                contribution = Pair{OccupationMonomial{S},D}[
                    transformed_monomial => transformed_coefficient
                ]
                _push_kernel_polynomial!(
                    out, transformed_sector, OccupationPolynomial{D,S}(contribution)
                )
            end
        end
    end

    return LoopQuotientedExpression{D,S,O,G,Ctx}(
        out, target_family(expression), parameters(expression), wigner_context(expression)
    )
end

"""
    quotient_loop_momenta(expression, workspace::LoopMomentumQuotientWorkspace)

Apply the exact physical dummy-loop quotient using caller-owned reusable scratch storage.

Construct `workspace` from an expression with the required graph/loop capacity and reuse it for
repeated quotient evaluations of that expression. The returned `LoopQuotientedExpression` is
identical to the one-shot `quotient_loop_momenta(expression)` result; only scratch ownership
differs.
"""
function quotient_loop_momenta(
    expression::OccupationReducedExpression{C,Boson,O,G,Ctx},
    workspace::LoopMomentumQuotientWorkspace,
) where {C<:Number,O,G,Ctx<:AbstractWignerContext}
    return _quotient_loop_momenta_graphcombinations(expression, workspace.storage)
end

function quotient_loop_momenta(
    expression::OccupationReducedExpression{C,Fermion,O,G,Ctx},
    workspace::LoopMomentumQuotientWorkspace,
) where {C<:Number,O,G,Ctx<:AbstractWignerContext}
    return _quotient_loop_momenta_graphcombinations(expression, workspace.storage)
end
