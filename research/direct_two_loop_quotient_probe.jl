using BenchmarkTools
using KeldyshContraction
using Test

import KeldyshContraction as KC

@qfields direct_two_loop_ϕ::Boson
@qfields direct_two_loop_ψ::Fermion

const TWO_LOOP_SUPPORT_KINDS = (:none, :shell, :pv, :mixed)
const TWO_LOOP_KINEMATIC_KINDS = (:none, :single, :pair, :mixed, :quartic)
const TWO_LOOP_OCCUPATION_KINDS = (:sparse, :repeated, :coupled, :dense)

function two_loop_support(
    ::Type{S}, family, k, q, r, kind::Symbol
) where {S<:KC.Statistics}
    energy(momentum) = KC.EnergyForm(KC.DispersionAtom(family, momentum))
    p = -k + q + r
    shells = KC.EnergyShell{S}[]
    pvs = KC.PrincipalValueSupport{S}[]
    if kind === :shell || kind === :mixed
        shell, _ = KC.energy_shell(energy(k) + energy(p) - energy(q) - energy(r))
        push!(shells, shell)
    end
    if kind === :pv || kind === :mixed
        pv, _ = KC.principal_value_support(energy(k) + energy(q) - energy(r))
        push!(pvs, pv)
    end
    return KC.FrequencySupport(shells, pvs)
end

function two_loop_kinematic(k, q, r, kind::Symbol)
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
        error("unknown two-loop kinematic fixture: $kind")
    end
    return KC.MomentumPolynomial(KC.MomentumMonomial(components), one(C))
end

function two_loop_monomial(
    ::Type{S}, family, k, q, r, kind::Symbol
) where {S<:KC.Statistics}
    momenta = if kind === :sparse
        [q]
    elseif kind === :repeated
        [q, q, r]
    elseif kind === :coupled
        [q + r, q - r]
    elseif kind === :dense
        [q, r, -k + q + r]
    else
        error("unknown two-loop occupation fixture: $kind")
    end
    return KC.OccupationMonomial(
        KC.OccupationAtom{S}[KC.OccupationAtom(family, momentum) for momentum in momenta]
    )
end

function two_loop_expression(
    ::Type{S}, family, support_kind::Symbol, kinematic_kind::Symbol, occupation_kind::Symbol
) where {S<:KC.Statistics}
    C = KC.ComplexRationals
    basis = KC.MomentumBasis(3)
    external = basis[1]
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    parameter = KC.ParameterMonomial(:direct_two_loop)
    sector = KC.ReducedCollisionSector{S}(
        parameter,
        basis,
        external,
        two_loop_kinematic(k, q, r, kinematic_kind),
        two_loop_support(S, family, k, q, r, support_kind),
    )
    monomial = two_loop_monomial(S, family, k, q, r, occupation_kind)
    polynomial = KC.OccupationPolynomial{C,S}([monomial => one(C)])
    return KC.OccupationReducedExpression{C,S,2,0,KC.HomogeneousWignerContext}(
        Dict(sector => polynomial), family, parameter, KC.HomogeneousWignerContext()
    )
end

@inline function direct_complex_isless(a, b)
    real(a) == real(b) || return real(a) < real(b)
    return imag(a) < imag(b)
end

function direct_momentum_polynomial_isless(a::KC.MomentumPolynomial, b::KC.MomentumPolynomial)
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        am, ac = a.terms[i]
        bm, bc = b.terms[i]
        am == bm || return isless(am, bm)
        ac == bc || return direct_complex_isless(ac, bc)
    end
    return length(a) < length(b)
end

function direct_energy_vector_isless(a, b)
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        ae = a[i].energy
        be = b[i].energy
        ae == be || return isless(ae, be)
    end
    return length(a) < length(b)
end

function direct_support_isless(a::KC.FrequencySupport, b::KC.FrequencySupport)
    if a.shells != b.shells
        return direct_energy_vector_isless(a.shells, b.shells)
    end
    return direct_energy_vector_isless(a.principal_values, b.principal_values)
end

function direct_candidate_isless(
    a_sector::KC.CollisionKernelSector{S},
    a_monomial::KC.OccupationMonomial{S},
    b_sector::KC.CollisionKernelSector{S},
    b_monomial::KC.OccupationMonomial{S},
) where {S<:KC.Statistics}
    a_monomial == b_monomial || return isless(a_monomial, b_monomial)
    ak = KC.kinematic_factor(a_sector)
    bk = KC.kinematic_factor(b_sector)
    ak == bk || return direct_momentum_polynomial_isless(ak, bk)
    as = KC.frequency_support(a_sector)
    bs = KC.frequency_support(b_sector)
    as == bs || return direct_support_isless(as, bs)
    return false
end

function two_loop_transforms(sector::KC.ReducedCollisionSector)
    basis = KC.momentum_basis(sector)
    length(basis) == 3 || error("direct two-loop quotient requires exactly two dummy loops")
    external = KC.external_wigner_momentum(sector)
    transforms = KC.LoopMomentumTransform[]
    shifts = [0, 0]
    for permutation in ([1, 2], [2, 1])
        for sign_mask in 0:3
            signs = [iszero(sign_mask & 1) ? 1 : -1, iszero(sign_mask & 2) ? 1 : -1]
            push!(
                transforms,
                KC.loop_permutation_transform(
                    basis, external, permutation, signs, shifts
                ),
            )
        end
    end
    return transforms
end

function direct_two_loop_quotient(
    expression::KC.OccupationReducedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:KC.Statistics,O,G,Ctx<:KC.AbstractWignerContext}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{D,S}}()

    for (sector, polynomial) in KC.occupation_reduced_terms(expression)
        transforms = two_loop_transforms(sector)
        for (occupation_monomial, occupation_coefficient) in polynomial
            for (kinematic_monomial, kinematic_coefficient) in KC.kinematic_factor(sector)
                atom_sector = KC._kinematic_atom_sector(sector, kinematic_monomial)
                best_sector = nothing
                best_monomial = nothing
                best_support_factor = one(Rational{Int})

                for transform in transforms
                    candidate_sector, support_factor = KC._transform_kernel_sector(
                        atom_sector, transform
                    )
                    candidate_monomial = KC.transform_loop_momenta(
                        occupation_monomial, transform
                    )
                    if best_sector === nothing || direct_candidate_isless(
                        candidate_sector,
                        candidate_monomial,
                        best_sector,
                        best_monomial,
                    )
                        best_sector = candidate_sector
                        best_monomial = candidate_monomial
                        best_support_factor = support_factor
                    end
                end

                coefficient =
                    convert(D, occupation_coefficient) *
                    convert(D, kinematic_coefficient) *
                    convert(D, best_support_factor)
                KC._push_kernel_polynomial!(
                    out,
                    best_sector,
                    KC.OccupationPolynomial{D,S}([best_monomial => coefficient]),
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

function nauty_two_loop_quotient(
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
                coefficient =
                    convert(D, occupation_coefficient) *
                    convert(D, kinematic_coefficient) *
                    convert(D, support_factor)
                KC._push_kernel_polynomial!(
                    out,
                    transformed_sector,
                    KC.OccupationPolynomial{D,S}([transformed_monomial => coefficient]),
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

function reflected_two_loop_expression(
    expression::KC.OccupationReducedExpression{C,S,O,G,Ctx},
    transform::KC.LoopMomentumTransform,
) where {C<:Number,S<:KC.Statistics,O,G,Ctx<:KC.AbstractWignerContext}
    sector, polynomial = only(KC.occupation_reduced_terms(expression))
    support, support_factor = KC._transform_frequency_support(
        KC.frequency_support(sector), transform
    )
    transformed_sector = KC.ReducedCollisionSector{S}(
        KC.parameters(sector),
        KC.momentum_basis(sector),
        KC.external_wigner_momentum(sector),
        KC.transform_loop_momenta(KC.kinematic_factor(sector), transform),
        support,
    )
    terms = Pair{KC.OccupationMonomial{S},C}[]
    for (monomial, coefficient) in polynomial
        push!(
            terms,
            KC.transform_loop_momenta(monomial, transform) =>
                coefficient * convert(C, support_factor),
        )
    end
    transformed_polynomial = KC.OccupationPolynomial{C,S}(terms)
    return KC.OccupationReducedExpression{C,S,O,G,Ctx}(
        Dict(transformed_sector => transformed_polynomial),
        KC.target_family(expression),
        KC.parameters(expression),
        KC.wigner_context(expression),
    )
end

function reduced_two_loop_sector(
    sector::KC.CollisionKernelSector{S}
) where {S<:KC.Statistics}
    return KC.ReducedCollisionSector{S}(
        KC.parameters(sector),
        KC.momentum_basis(sector),
        KC.external_wigner_momentum(sector),
        KC.kinematic_factor(sector),
        KC.frequency_support(sector),
    )
end

function two_loop_gc_workspace(
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
    for (kinematic_monomial, _) in KC.kinematic_factor(sector)
        graph_capacity += 1
        for component in kinematic_monomial
            graph_capacity += KC._projective_component_vertex_count(component)
        end
    end
    return KC._LoopGCQuotientWorkspace(graph_capacity, nloops)
end

function gc_requotient_two_loop_terms(
    expression::KC.LoopQuotientedExpression{C,S}
) where {C<:Number,S<:KC.Statistics}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{D,S}}()
    for (sector, polynomial) in KC.loop_quotient_terms(expression)
        source_sector = reduced_two_loop_sector(sector)
        for (monomial, coefficient) in polynomial
            workspace = two_loop_gc_workspace(source_sector, monomial)
            transform = KC._graphcombinations_projective_canonical_loop_transform(
                source_sector, monomial, workspace
            )
            transformed_sector, support_factor = KC._transform_kernel_sector(
                source_sector, transform
            )
            transformed_monomial = KC.transform_loop_momenta(monomial, transform)
            coefficient2 = convert(D, coefficient) * convert(D, support_factor)
            KC._push_kernel_polynomial!(
                out,
                transformed_sector,
                KC.OccupationPolynomial{D,S}([transformed_monomial => coefficient2]),
            )
        end
    end
    return out
end

function two_loop_cases()
    out = Pair{String,Any}[]
    for (S, family) in ((Boson, direct_two_loop_ϕ), (Fermion, direct_two_loop_ψ))
        statistics = S === Boson ? "boson" : "fermion"
        for support in TWO_LOOP_SUPPORT_KINDS,
            kinematic in TWO_LOOP_KINEMATIC_KINDS,
            occupation in TWO_LOOP_OCCUPATION_KINDS
            label = "$statistics/$support/$kinematic/$occupation"
            push!(
                out,
                label => two_loop_expression(S, family, support, kinematic, occupation),
            )
        end
    end
    return out
end

all_two_loop_cases = two_loop_cases()
two_loop_semantic_failures = String[]
two_loop_invariance_failures = String[]
two_loop_time_regressions = Tuple{String,Float64}[]
two_loop_memory_regressions = Tuple{String,Int,Int}[]
worst_direct_nauty_ratio = ("", 0.0)
worst_direct_gc_ratio = ("", 0.0)

@testset "direct two-loop semantics" begin
    for (label, expression) in all_two_loop_cases
        direct = direct_two_loop_quotient(expression)
        nauty = nauty_two_loop_quotient(expression)
        gc = KC.quotient_loop_momenta(expression)
        oracle = KC.loop_quotient_terms(gc)
        exact =
            gc_requotient_two_loop_terms(direct) == oracle &&
            gc_requotient_two_loop_terms(nauty) == oracle
        exact || push!(two_loop_semantic_failures, label)
        @test exact

        sector = first(first(KC.occupation_reduced_terms(expression)))
        reference = KC.loop_quotient_terms(direct)
        invariant = true
        for transform in two_loop_transforms(sector)
            transformed = reflected_two_loop_expression(expression, transform)
            if KC.loop_quotient_terms(direct_two_loop_quotient(transformed)) != reference
                invariant = false
                break
            end
        end
        invariant || push!(two_loop_invariance_failures, label)
        @test invariant
    end
end

@testset "direct two-loop performance" begin
    global worst_direct_nauty_ratio, worst_direct_gc_ratio
    for (label, expression) in all_two_loop_cases
        direct_two_loop_quotient(expression)
        nauty_two_loop_quotient(expression)
        KC.quotient_loop_momenta(expression)

        direct_trial = @benchmark direct_two_loop_quotient($expression) samples = 21 evals = 1
        nauty_trial = @benchmark nauty_two_loop_quotient($expression) samples = 21 evals = 1
        gc_trial = @benchmark KC.quotient_loop_momenta($expression) samples = 21 evals = 1
        direct = median(direct_trial)
        nauty = median(nauty_trial)
        gc = median(gc_trial)

        if direct.time > nauty.time || direct.memory > nauty.memory
            direct_trial = @benchmark direct_two_loop_quotient($expression) samples = 61 evals = 1
            nauty_trial = @benchmark nauty_two_loop_quotient($expression) samples = 61 evals = 1
            direct = median(direct_trial)
            nauty = median(nauty_trial)
        end

        direct_nauty_ratio = direct.time / nauty.time
        direct_gc_ratio = direct.time / gc.time
        direct_nauty_ratio > worst_direct_nauty_ratio[2] &&
            (worst_direct_nauty_ratio = (label, direct_nauty_ratio))
        direct_gc_ratio > worst_direct_gc_ratio[2] &&
            (worst_direct_gc_ratio = (label, direct_gc_ratio))
        direct.time > nauty.time &&
            push!(two_loop_time_regressions, (label, direct_nauty_ratio))
        direct.memory > nauty.memory &&
            push!(two_loop_memory_regressions, (label, direct.memory, nauty.memory))
    end

    println("semantic failures: ", length(two_loop_semantic_failures))
    println("signed-permutation invariance failures: ", length(two_loop_invariance_failures))
    println("direct two-loop time regressions vs Nauty: ", length(two_loop_time_regressions))
    println(
        "direct two-loop memory regressions vs Nauty: ",
        length(two_loop_memory_regressions),
    )
    println(
        "worst direct/Nauty ratio: ",
        round(worst_direct_nauty_ratio[2]; digits=3),
        "x at ",
        worst_direct_nauty_ratio[1],
    )
    println(
        "worst direct/GC ratio: ",
        round(worst_direct_gc_ratio[2]; digits=3),
        "x at ",
        worst_direct_gc_ratio[1],
    )
    @test isempty(two_loop_time_regressions)
    @test isempty(two_loop_memory_regressions)
end
