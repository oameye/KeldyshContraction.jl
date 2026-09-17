using BenchmarkTools
using KeldyshContraction
using Test

import KeldyshContraction as KC

@qfields gc_corpus_ϕ::Boson
@qfields gc_corpus_ψ::Fermion

const CORPUS_LOOP_COUNTS = 1:5
const CORPUS_SUPPORTS = (:none, :shell, :pv, :mixed)
const CORPUS_KINEMATICS = (:none, :single, :pair, :mixed, :quartic)
const CORPUS_OCCUPATIONS = (:sparse, :repeated, :coupled, :dense)

function _corpus_support(
    ::Type{S}, family, k, q, r, p, kind::Symbol
) where {S<:KC.Statistics}
    energy(momentum) = KC.EnergyForm(KC.DispersionAtom(family, momentum))
    shells = KC.EnergyShell{S}[]
    principal_values = KC.PrincipalValueSupport{S}[]

    if kind === :shell || kind === :mixed
        shell, _ = KC.energy_shell(energy(k) + energy(p) - energy(q) - energy(r))
        push!(shells, shell)
    end
    if kind === :pv || kind === :mixed
        pv, _ = KC.principal_value_support(energy(k) + energy(q) - energy(r))
        push!(principal_values, pv)
    end
    return KC.FrequencySupport(shells, principal_values)
end

function _corpus_kinematic(k, q, r, kind::Symbol)
    C = KC.ComplexRationals
    components = if kind === :none
        KC.MomentumComponent[]
    elseif kind === :single
        [KC.MomentumComponent(q, :x)]
    elseif kind === :pair
        [KC.MomentumComponent(q, :x), KC.MomentumComponent(r, :x)]
    elseif kind === :mixed
        [KC.MomentumComponent(k + q - r, :x), KC.MomentumComponent(q + r, :y)]
    elseif kind === :quartic
        qx = KC.MomentumComponent(q, :x)
        rx = KC.MomentumComponent(r, :x)
        [qx, qx, rx, rx]
    else
        error("unknown kinematic corpus kind: $kind")
    end
    return KC.MomentumPolynomial(KC.MomentumMonomial(components), one(C))
end

function _corpus_monomial(
    ::Type{S}, family, k, loops, q, r, kind::Symbol
) where {S<:KC.Statistics}
    momenta = if kind === :sparse
        [q]
    elseif kind === :repeated
        [q, q, r]
    elseif kind === :coupled
        [q + r, q - r]
    elseif kind === :dense
        total = reduce(+, loops)
        vcat(loops, [-k + total])
    else
        error("unknown occupation corpus kind: $kind")
    end
    atoms = KC.OccupationAtom{S}[
        KC.OccupationAtom(family, momentum) for momentum in momenta
    ]
    return KC.OccupationMonomial(atoms)
end

function corpus_expression(
    ::Type{S},
    family,
    nloops::Int,
    support_kind::Symbol,
    kinematic_kind::Symbol,
    occupation_kind::Symbol,
) where {S<:KC.Statistics}
    C = KC.ComplexRationals
    basis = KC.MomentumBasis(nloops + 1)
    external = basis[1]
    k = KC.basis_momentum(basis, 1)
    loops = [KC.basis_momentum(basis, i) for i in 2:(nloops + 1)]
    q = loops[1]
    r = loops[min(2, nloops)]
    p = -k + q + r
    parameter = KC.ParameterMonomial(:gc_corpus)

    support = _corpus_support(S, family, k, q, r, p, support_kind)
    kinematic = _corpus_kinematic(k, q, r, kinematic_kind)
    sector = KC.ReducedCollisionSector{S}(parameter, basis, external, kinematic, support)
    monomial = _corpus_monomial(S, family, k, loops, q, r, occupation_kind)
    polynomial = KC.OccupationPolynomial{C,S}([monomial => one(C)])

    return KC.OccupationReducedExpression{C,S,2,0,KC.HomogeneousWignerContext}(
        Dict(sector => polynomial), family, parameter, KC.HomogeneousWignerContext()
    )
end

function nauty_quotient_loop_momenta(
    expression::KC.OccupationReducedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:KC.Statistics,O,G,Ctx<:KC.AbstractWignerContext}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{D,S}}()

    for (sector, polynomial) in KC.occupation_reduced_terms(expression)
        for (occupation_monomial, occupation_coefficient) in polynomial
            for (kinematic_monomial, kinematic_coefficient) in KC.kinematic_factor(sector)
                atom_sector = KC._kinematic_atom_sector(sector, kinematic_monomial)
                transform = KC._projective_canonical_loop_transform(
                    atom_sector, occupation_monomial
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

function _reduced_sector(sector::KC.CollisionKernelSector{S}) where {S<:KC.Statistics}
    return KC.ReducedCollisionSector{S}(
        KC.parameters(sector),
        KC.momentum_basis(sector),
        KC.external_wigner_momentum(sector),
        KC.kinematic_factor(sector),
        KC.frequency_support(sector),
    )
end

function _gc_oracle_workspace(
    sector::KC.ReducedCollisionSector, monomial::KC.OccupationMonomial
)
    basis = KC.momentum_basis(sector)
    external_index = KC._external_basis_index_noalloc(
        basis, KC.external_wigner_momentum(sector)
    )
    nloops = length(basis) - 1
    graph_capacity =
        3 + 3 * nloops + KC._projective_support_vertex_count(sector, external_index)

    for atom in monomial
        graph_capacity += 1 + KC._loop_incidence_vertex_count(atom.momentum, external_index)
    end

    kinematic_capacity = 0
    for (kinematic_monomial, _) in KC.kinematic_factor(sector)
        vertices = 1
        for component in kinematic_monomial
            vertices += KC._projective_component_vertex_count(component)
        end
        kinematic_capacity = max(kinematic_capacity, vertices)
    end
    graph_capacity += kinematic_capacity
    return KC._LoopGCQuotientWorkspace(graph_capacity, nloops)
end

function gc_requotient_terms(
    expression::KC.LoopQuotientedExpression{C,S}
) where {C<:Number,S<:KC.Statistics}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{D,S}}()

    for (sector, polynomial) in KC.loop_quotient_terms(expression)
        source_sector = _reduced_sector(sector)
        for (monomial, coefficient) in polynomial
            workspace = _gc_oracle_workspace(source_sector, monomial)
            transform = KC._graphcombinations_projective_canonical_loop_transform(
                source_sector, monomial, workspace
            )
            transformed_sector, support_factor = KC._transform_kernel_sector(
                source_sector, transform
            )
            transformed_monomial = KC.transform_loop_momenta(monomial, transform)
            transformed_coefficient = convert(D, coefficient) * convert(D, support_factor)
            contribution = Pair{KC.OccupationMonomial{S},D}[transformed_monomial => transformed_coefficient]
            KC._push_kernel_polynomial!(
                out, transformed_sector, KC.OccupationPolynomial{D,S}(contribution)
            )
        end
    end
    return out
end

function _case_label(S, nloops, support, kinematic, occupation)
    statistics = S === Boson ? "boson" : "fermion"
    return "$statistics/$nloops-loop/$support/$kinematic/$occupation"
end

function corpus_cases()
    cases = Pair{String,Any}[]
    for (S, family) in ((Boson, gc_corpus_ϕ), (Fermion, gc_corpus_ψ))
        for nloops in CORPUS_LOOP_COUNTS
            for support in CORPUS_SUPPORTS
                for kinematic in CORPUS_KINEMATICS
                    for occupation in CORPUS_OCCUPATIONS
                        label = _case_label(S, nloops, support, kinematic, occupation)
                        expression = corpus_expression(
                            S, family, nloops, support, kinematic, occupation
                        )
                        push!(cases, label => expression)
                    end
                end
            end
        end
    end
    return cases
end

function _measure_pair_15(expression)
    nauty_quotient_loop_momenta(expression)
    KC.quotient_loop_momenta(expression)
    nauty_trial = @benchmark nauty_quotient_loop_momenta($expression) samples = 15 evals = 1
    gc_trial = @benchmark KC.quotient_loop_momenta($expression) samples = 15 evals = 1
    return median(nauty_trial), median(gc_trial)
end

function _measure_pair_41(expression)
    nauty_trial = @benchmark nauty_quotient_loop_momenta($expression) samples = 41 evals = 1
    gc_trial = @benchmark KC.quotient_loop_momenta($expression) samples = 41 evals = 1
    return median(nauty_trial), median(gc_trial)
end

function _confirmed_measurement(expression)
    nauty, gc = _measure_pair_15(expression)
    if gc.time > nauty.time || gc.memory > nauty.memory
        nauty, gc = _measure_pair_41(expression)
    end
    return nauty, gc
end

cases = corpus_cases()
println("GC/Nauty loop corpus: ", length(cases), " generated production-shaped cases")

semantic_failures = String[]
time_regressions = Tuple{String,Float64}[]
memory_regressions = Tuple{String,Int,Int}[]
worst_ratio = ("", 0.0)
best_ratio = ("", Inf)

@testset "GC/Nauty loop corpus semantics" begin
    for (label, expression) in cases
        nauty = nauty_quotient_loop_momenta(expression)
        gc = KC.quotient_loop_momenta(expression)
        gc_terms = KC.loop_quotient_terms(gc)
        exact =
            gc_requotient_terms(gc) == gc_terms && gc_requotient_terms(nauty) == gc_terms
        if !exact
            push!(semantic_failures, label)
        end
        @test exact
    end
end

@testset "GC/Nauty loop corpus performance" begin
    for (index, (label, expression)) in enumerate(cases)
        nauty, gc = _confirmed_measurement(expression)
        ratio = gc.time / nauty.time
        if ratio > worst_ratio[2]
            worst_ratio = (label, ratio)
        end
        if ratio < best_ratio[2]
            best_ratio = (label, ratio)
        end
        if gc.time > nauty.time
            push!(time_regressions, (label, ratio))
        end
        if gc.memory > nauty.memory
            push!(memory_regressions, (label, gc.memory, nauty.memory))
        end
        println(
            lpad(index, 3),
            "/",
            length(cases),
            " ",
            label,
            ": time=",
            round(gc.time / 1.0e3; digits=2),
            "/",
            round(nauty.time / 1.0e3; digits=2),
            " μs GC/Nauty (",
            round(ratio; digits=3),
            "x), memory=",
            gc.memory,
            "/",
            nauty.memory,
            " B",
        )
    end

    println("semantic failures: ", length(semantic_failures))
    println("confirmed time regressions: ", length(time_regressions))
    println("memory regressions: ", length(memory_regressions))
    println("best time ratio: ", round(best_ratio[2]; digits=3), "x at ", best_ratio[1])
    println("worst time ratio: ", round(worst_ratio[2]; digits=3), "x at ", worst_ratio[1])

    if !isempty(time_regressions)
        println("confirmed GC time regressions:")
        for (label, ratio) in sort(time_regressions; by=last, rev=true)
            println("  ", label, ": ", round(ratio; digits=3), "x")
        end
    end
    if !isempty(memory_regressions)
        println("GC memory regressions:")
        for (label, gc_memory, nauty_memory) in memory_regressions
            println("  ", label, ": ", gc_memory, " B vs ", nauty_memory, " B")
        end
    end

    @test isempty(time_regressions)
    @test isempty(memory_regressions)
end
