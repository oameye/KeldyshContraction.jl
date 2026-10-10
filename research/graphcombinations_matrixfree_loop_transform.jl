using BenchmarkTools
using Test

import GraphCombinations as GC
import KeldyshContraction as KC

include(ENV["KC_CORPUS_DEFS"])

const RESIDUAL_TWO_LOOP_LABELS = Set([
    "boson/2-loop/none/none/sparse",
    "boson/2-loop/none/none/dense",
    "boson/2-loop/none/single/sparse",
    "boson/2-loop/none/pair/dense",
    "boson/2-loop/shell/none/dense",
    "boson/2-loop/shell/pair/dense",
    "boson/2-loop/pv/none/sparse",
    "fermion/2-loop/none/none/sparse",
    "fermion/2-loop/none/none/dense",
    "fermion/2-loop/shell/none/dense",
    "fermion/2-loop/shell/pair/dense",
    "fermion/2-loop/pv/none/sparse",
    "fermion/2-loop/pv/single/sparse",
])

function _prepare_matrixfree_gauge!(
    sector::KC.ReducedCollisionSector{S},
    monomial::KC.OccupationMonomial{S},
    workspace::KC._LoopGCQuotientWorkspace,
) where {S<:KC.Statistics}
    basis = KC.momentum_basis(sector)
    external = KC.external_wigner_momentum(sector)
    nloops = length(basis) - 1
    KC._prepare_loop_gc_workspace!(workspace, nloops)
    external_index = KC._write_loop_basis_indices!(workspace, basis, external)

    builder = workspace.builder
    root = KC._add_loop_vertex!(builder, KC._loop_graph_color(1))
    for slot in 1:nloops
        pair = KC._add_loop_vertex!(builder, KC._loop_graph_color(2))
        positive = KC._add_loop_vertex!(builder, KC._loop_graph_color(3))
        negative = KC._add_loop_vertex!(builder, KC._loop_graph_color(3))
        workspace.pair_vertices[slot] = pair
        workspace.positive_vertices[slot] = positive
        workspace.negative_vertices[slot] = negative
        KC._add_loop_edge!(builder, root, pair)
        KC._add_loop_edge!(builder, pair, positive)
        KC._add_loop_edge!(builder, pair, negative)
    end

    KC._add_projective_loop_semantics!(
        builder,
        root,
        sector,
        monomial,
        external_index,
        workspace.loop_indices,
        workspace.positive_vertices,
        workspace.negative_vertices,
        workspace.loop_incidences,
    )

    buffer = KC._canonicalize_loop_builder!(workspace.canonicalization, builder)
    KC._order_loop_slots!(workspace, buffer, nloops)
    @inbounds for new_slot in 1:nloops
        old_slot = workspace.ordered_slots[new_slot]
        workspace.loop_permutation[old_slot] = new_slot
        positive_rank = GC.canonical_rank(buffer, workspace.positive_vertices[old_slot])
        negative_rank = GC.canonical_rank(buffer, workspace.negative_vertices[old_slot])
        workspace.loop_signs[old_slot] = positive_rank < negative_rank ? 1 : -1
    end
    return external_index
end

function _matrixfree_momentum(
    momentum::KC.LinearMomentum, workspace::KC._LoopGCQuotientWorkspace, external_index::Int
)
    n = length(momentum)
    n == length(workspace.loop_indices) + 1 || throw(
        DimensionMismatch("momentum and signed loop permutation use different basis sizes"),
    )
    coefficients = Vector{KC.MomentumCoefficient}(undef, n)
    coefficients[external_index] = momentum[external_index]
    @inbounds for old_slot in eachindex(workspace.loop_permutation)
        old_index = workspace.loop_indices[old_slot]
        new_slot = workspace.loop_permutation[old_slot]
        new_index = workspace.loop_indices[new_slot]
        coefficients[new_index] = workspace.loop_signs[old_slot] * momentum[old_index]
    end
    return KC.LinearMomentum(coefficients)
end

function _matrixfree_occupation_monomial(
    monomial::KC.OccupationMonomial{S},
    workspace::KC._LoopGCQuotientWorkspace,
    external_index::Int,
) where {S<:KC.Statistics}
    atoms = KC.OccupationAtom{S}[
        KC.OccupationAtom{S}(
            atom.family, _matrixfree_momentum(atom.momentum, workspace, external_index)
        ) for atom in monomial
    ]
    return KC.OccupationMonomial(atoms)
end

function _matrixfree_momentum_monomial(
    monomial::KC.MomentumMonomial,
    workspace::KC._LoopGCQuotientWorkspace,
    external_index::Int,
)
    factors = KC.MomentumComponent[
        KC.MomentumComponent(
            _matrixfree_momentum(component.momentum, workspace, external_index),
            component.axis,
        ) for component in monomial
    ]
    return KC.MomentumMonomial(factors)
end

function _matrixfree_momentum_polynomial(
    polynomial::KC.MomentumPolynomial{C},
    workspace::KC._LoopGCQuotientWorkspace,
    external_index::Int,
) where {C<:Number}
    terms = Pair{KC.MomentumMonomial,C}[]
    sizehint!(terms, length(polynomial))
    for (monomial, coefficient) in polynomial
        transformed = _matrixfree_momentum_monomial(monomial, workspace, external_index)
        normalized, factor, nonzero = KC._projective_kinematic_monomial(transformed)
        nonzero || continue
        push!(terms, normalized => coefficient * convert(C, factor))
    end
    return KC.MomentumPolynomial{C}(terms)
end

function _matrixfree_energy_form(
    form::KC.EnergyForm{S}, workspace::KC._LoopGCQuotientWorkspace, external_index::Int
) where {S<:KC.Statistics}
    KC.energy_basis_size(form) == length(workspace.loop_indices) + 1 || throw(
        DimensionMismatch(
            "energy form and signed loop permutation use different basis sizes"
        ),
    )
    terms = Pair{KC.DispersionAtom{S},KC.EnergyCoefficient}[
        KC.DispersionAtom{S}(
            atom.family, _matrixfree_momentum(atom.momentum, workspace, external_index)
        ) => coefficient for (atom, coefficient) in form
    ]
    return KC.EnergyForm(KC.energy_basis_size(form), terms)
end

function _matrixfree_frequency_support(
    support::KC.FrequencySupport{S},
    workspace::KC._LoopGCQuotientWorkspace,
    external_index::Int,
) where {S<:KC.Statistics}
    shells = KC.EnergyShell{S}[]
    principal_values = KC.PrincipalValueSupport{S}[]
    factor = one(KC.MomentumCoefficient)
    sizehint!(shells, length(support.shells))
    sizehint!(principal_values, length(support.principal_values))

    for shell in support.shells
        transformed, shell_factor = KC.energy_shell(
            _matrixfree_energy_form(shell.energy, workspace, external_index)
        )
        push!(shells, transformed)
        factor *= shell_factor
    end
    for principal_value in support.principal_values
        transformed, pv_factor = KC.principal_value_support(
            _matrixfree_energy_form(principal_value.energy, workspace, external_index)
        )
        push!(principal_values, transformed)
        factor *= pv_factor
    end
    return KC.FrequencySupport(shells, principal_values), factor
end

function _matrixfree_sector(
    sector::KC.ReducedCollisionSector{S},
    workspace::KC._LoopGCQuotientWorkspace,
    external_index::Int,
) where {S<:KC.Statistics}
    support, factor = _matrixfree_frequency_support(
        KC.frequency_support(sector), workspace, external_index
    )
    transformed = KC.CollisionKernelSector{S}(
        KC.parameters(sector),
        KC.momentum_basis(sector),
        KC.external_wigner_momentum(sector),
        _matrixfree_momentum_polynomial(
            KC.kinematic_factor(sector), workspace, external_index
        ),
        support,
    )
    return transformed, factor
end

function prepared_dense_gc_quotient(
    expression::KC.OccupationReducedExpression{C,S,O,G,Ctx},
    workspace::KC._LoopGCQuotientWorkspace,
) where {C<:Number,S<:KC.Statistics,O,G,Ctx<:KC.AbstractWignerContext}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{D,S}}()
    for (sector, polynomial) in KC.occupation_reduced_terms(expression)
        for (occupation_monomial, occupation_coefficient) in polynomial
            for (kinematic_monomial, kinematic_coefficient) in KC.kinematic_factor(sector)
                atom_sector = KC._kinematic_atom_sector(sector, kinematic_monomial)
                transform = KC._graphcombinations_projective_canonical_loop_transform(
                    atom_sector, occupation_monomial, workspace
                )
                transformed_sector, support_factor = KC._transform_kernel_sector(
                    atom_sector, transform
                )
                transformed_monomial = KC.transform_loop_momenta(
                    occupation_monomial, transform
                )
                transformed_coefficient =
                    convert(D, occupation_coefficient) *
                    convert(D, kinematic_coefficient) *
                    convert(D, support_factor)
                contribution = Pair{KC.OccupationMonomial{S},D}[transformed_monomial => transformed_coefficient]
                KC._push_kernel_polynomial!(
                    out, transformed_sector, KC.OccupationPolynomial{D,S}(contribution)
                )
            end
        end
    end
    return KC.LoopQuotientedExpression{D,S,O,G,Ctx}(
        out,
        KC.target_family(expression),
        KC.parameters(expression),
        KC.wigner_context(expression),
    )
end

function matrixfree_gc_quotient(
    expression::KC.OccupationReducedExpression{C,S,O,G,Ctx},
    workspace::KC._LoopGCQuotientWorkspace,
) where {C<:Number,S<:KC.Statistics,O,G,Ctx<:KC.AbstractWignerContext}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{D,S}}()
    for (sector, polynomial) in KC.occupation_reduced_terms(expression)
        for (occupation_monomial, occupation_coefficient) in polynomial
            for (kinematic_monomial, kinematic_coefficient) in KC.kinematic_factor(sector)
                atom_sector = KC._kinematic_atom_sector(sector, kinematic_monomial)
                external_index = _prepare_matrixfree_gauge!(
                    atom_sector, occupation_monomial, workspace
                )
                transformed_sector, support_factor = _matrixfree_sector(
                    atom_sector, workspace, external_index
                )
                transformed_monomial = _matrixfree_occupation_monomial(
                    occupation_monomial, workspace, external_index
                )
                transformed_coefficient =
                    convert(D, occupation_coefficient) *
                    convert(D, kinematic_coefficient) *
                    convert(D, support_factor)
                contribution = Pair{KC.OccupationMonomial{S},D}[transformed_monomial => transformed_coefficient]
                KC._push_kernel_polynomial!(
                    out, transformed_sector, KC.OccupationPolynomial{D,S}(contribution)
                )
            end
        end
    end
    return KC.LoopQuotientedExpression{D,S,O,G,Ctx}(
        out,
        KC.target_family(expression),
        KC.parameters(expression),
        KC.wigner_context(expression),
    )
end

function _workspace(expression)
    graph_capacity, loop_capacity = KC._loop_gc_capacities(expression)
    return KC._LoopGCQuotientWorkspace(graph_capacity, loop_capacity)
end

function residual_cases()
    return [
        label => expression for
        (label, expression) in corpus_cases() if label in RESIDUAL_TWO_LOOP_LABELS
    ]
end

cases = residual_cases()
workspaces = [_workspace(expression) for (_, expression) in cases]
println("matrix-free packed-GC residual corpus: ", length(cases), " cases")

@testset "matrix-free signed-permutation semantics" begin
    @test length(cases) == length(RESIDUAL_TWO_LOOP_LABELS)
    for (index, (_, expression)) in enumerate(cases)
        workspace = workspaces[index]
        @test KC.loop_quotient_terms(matrixfree_gc_quotient(expression, workspace)) ==
            KC.loop_quotient_terms(prepared_dense_gc_quotient(expression, workspace))
    end
end

time_regressions = Tuple{String,Float64}[]
memory_regressions = Tuple{String,Int,Int}[]
dense_ratios = Tuple{String,Float64}[]

@testset "matrix-free signed-permutation residual performance" begin
    for (index, (label, expression)) in enumerate(cases)
        workspace = workspaces[index]
        nauty_quotient_loop_momenta(expression)
        prepared_dense_gc_quotient(expression, workspace)
        matrixfree_gc_quotient(expression, workspace)

        nauty_trial = @benchmark nauty_quotient_loop_momenta($expression) samples = 41 evals =
            1
        dense_trial = @benchmark prepared_dense_gc_quotient($expression, $workspace) samples =
            41 evals = 1
        matrixfree_trial = @benchmark matrixfree_gc_quotient($expression, $workspace) samples =
            41 evals = 1
        nauty = median(nauty_trial)
        dense = median(dense_trial)
        matrixfree = median(matrixfree_trial)
        nauty_ratio = matrixfree.time / nauty.time
        dense_ratio = matrixfree.time / dense.time
        push!(dense_ratios, (label, dense_ratio))
        matrixfree.time > nauty.time && push!(time_regressions, (label, nauty_ratio))
        matrixfree.memory > nauty.memory &&
            push!(memory_regressions, (label, matrixfree.memory, nauty.memory))
        println(
            label,
            ": matrixfree=",
            round(matrixfree.time / 1.0e3; digits=2),
            " μs/",
            matrixfree.memory,
            " B; dense=",
            round(dense.time / 1.0e3; digits=2),
            " μs/",
            dense.memory,
            " B; Nauty=",
            round(nauty.time / 1.0e3; digits=2),
            " μs/",
            nauty.memory,
            " B; ratios=",
            round(dense_ratio; digits=3),
            "x dense, ",
            round(nauty_ratio; digits=3),
            "x Nauty",
        )
    end

    println("confirmed time regressions: ", length(time_regressions))
    println("memory regressions: ", length(memory_regressions))
    println("matrix-free/dense ratios:")
    for (label, ratio) in sort(dense_ratios; by=last, rev=true)
        println("  ", label, ": ", round(ratio; digits=3), "x")
    end
    if !isempty(time_regressions)
        println("remaining matrix-free/Nauty time regressions:")
        for (label, ratio) in sort(time_regressions; by=last, rev=true)
            println("  ", label, ": ", round(ratio; digits=3), "x")
        end
    end
    if !isempty(memory_regressions)
        println("remaining matrix-free/Nauty memory regressions:")
        for (label, matrixfree_memory, nauty_memory) in memory_regressions
            println("  ", label, ": ", matrixfree_memory, " B vs ", nauty_memory, " B")
        end
    end

    @test isempty(time_regressions)
    @test isempty(memory_regressions)
end
