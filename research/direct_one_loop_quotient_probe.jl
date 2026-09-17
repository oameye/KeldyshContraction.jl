using BenchmarkTools
using KeldyshContraction
using Test

import KeldyshContraction as KC

@qfields direct_loop_ϕ::Boson
@qfields direct_loop_ψ::Fermion

const SUPPORT_KINDS = (:none, :shell, :pv, :mixed)
const KINEMATIC_KINDS = (:none, :single, :pair, :mixed, :quartic)
const OCCUPATION_KINDS = (:sparse, :repeated, :coupled, :dense)

function support_fixture(::Type{S}, family, k, q, kind::Symbol) where {S<:KC.Statistics}
    energy(momentum) = KC.EnergyForm(KC.DispersionAtom(family, momentum))
    shells = KC.EnergyShell{S}[]
    pvs = KC.PrincipalValueSupport{S}[]
    p = -k + 2q
    if kind === :shell || kind === :mixed
        shell, _ = KC.energy_shell(energy(k) + energy(p) - 2 * energy(q))
        push!(shells, shell)
    end
    if kind === :pv || kind === :mixed
        pv, _ = KC.principal_value_support(energy(k) + energy(q) - energy(p))
        push!(pvs, pv)
    end
    return KC.FrequencySupport(shells, pvs)
end

function kinematic_fixture(k, q, kind::Symbol)
    C = KC.ComplexRationals
    components = if kind === :none
        KC.MomentumComponent[]
    elseif kind === :single
        [KC.MomentumComponent(q, :x)]
    elseif kind === :pair
        [KC.MomentumComponent(q, :x), KC.MomentumComponent(q, :y)]
    elseif kind === :mixed
        [KC.MomentumComponent(k + q, :x), KC.MomentumComponent(k - q, :y)]
    elseif kind === :quartic
        qx = KC.MomentumComponent(q, :x)
        qy = KC.MomentumComponent(q, :y)
        [qx, qx, qy, qy]
    else
        error("unknown kinematic fixture: $kind")
    end
    return KC.MomentumPolynomial(KC.MomentumMonomial(components), one(C))
end

function occupation_fixture(::Type{S}, family, k, q, kind::Symbol) where {S<:KC.Statistics}
    momenta = if kind === :sparse
        [q]
    elseif kind === :repeated
        [q, q, k - q]
    elseif kind === :coupled
        [k + q, k - q]
    elseif kind === :dense
        [q, k - q, k + q, 2q - k]
    else
        error("unknown occupation fixture: $kind")
    end
    return KC.OccupationMonomial(
        KC.OccupationAtom{S}[KC.OccupationAtom(family, momentum) for momentum in momenta]
    )
end

function one_loop_expression(
    ::Type{S}, family, support_kind::Symbol, kinematic_kind::Symbol, occupation_kind::Symbol
) where {S<:KC.Statistics}
    C = KC.ComplexRationals
    basis = KC.MomentumBasis(2)
    external = basis[1]
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    parameter = KC.ParameterMonomial(:direct_loop)
    sector = KC.ReducedCollisionSector{S}(
        parameter,
        basis,
        external,
        kinematic_fixture(k, q, kinematic_kind),
        support_fixture(S, family, k, q, support_kind),
    )
    monomial = occupation_fixture(S, family, k, q, occupation_kind)
    polynomial = KC.OccupationPolynomial{C,S}([monomial => one(C)])
    return KC.OccupationReducedExpression{C,S,2,0,KC.HomogeneousWignerContext}(
        Dict(sector => polynomial), family, parameter, KC.HomogeneousWignerContext()
    )
end

function nauty_quotient(
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

@inline function complex_isless(a, b)
    real(a) == real(b) || return real(a) < real(b)
    return imag(a) < imag(b)
end

function momentum_polynomial_isless(a::KC.MomentumPolynomial, b::KC.MomentumPolynomial)
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        am, ac = a.terms[i]
        bm, bc = b.terms[i]
        am == bm || return isless(am, bm)
        ac == bc || return complex_isless(ac, bc)
    end
    return length(a) < length(b)
end

function energy_vector_isless(a, b)
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        ae = a[i].energy
        be = b[i].energy
        ae == be || return isless(ae, be)
    end
    return length(a) < length(b)
end

function support_isless(a::KC.FrequencySupport, b::KC.FrequencySupport)
    if a.shells != b.shells
        return energy_vector_isless(a.shells, b.shells)
    end
    return energy_vector_isless(a.principal_values, b.principal_values)
end

function candidate_isless(
    a_sector::KC.CollisionKernelSector{S},
    a_monomial::KC.OccupationMonomial{S},
    b_sector::KC.CollisionKernelSector{S},
    b_monomial::KC.OccupationMonomial{S},
) where {S<:KC.Statistics}
    a_monomial == b_monomial || return isless(a_monomial, b_monomial)
    ak = KC.kinematic_factor(a_sector)
    bk = KC.kinematic_factor(b_sector)
    ak == bk || return momentum_polynomial_isless(ak, bk)
    as = KC.frequency_support(a_sector)
    bs = KC.frequency_support(b_sector)
    as == bs || return support_isless(as, bs)
    return false
end

function one_loop_reflection(sector::KC.ReducedCollisionSector)
    basis = KC.momentum_basis(sector)
    length(basis) == 2 || error("direct one-loop quotient requires exactly one dummy loop")
    external = KC.external_wigner_momentum(sector)
    return KC.loop_permutation_transform(basis, external, [1], [-1], [0])
end

function untransformed_sector(sector::KC.ReducedCollisionSector{S}) where {S<:KC.Statistics}
    return KC.CollisionKernelSector{S}(
        KC.parameters(sector),
        KC.momentum_basis(sector),
        KC.external_wigner_momentum(sector),
        KC.kinematic_factor(sector),
        KC.frequency_support(sector),
    )
end

function direct_one_loop_quotient(
    expression::KC.OccupationReducedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:KC.Statistics,O,G,Ctx<:KC.AbstractWignerContext}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{D,S}}()
    for (sector, polynomial) in KC.occupation_reduced_terms(expression)
        reflection = one_loop_reflection(sector)
        for (occupation_monomial, occupation_coefficient) in polynomial
            for (kinematic_monomial, kinematic_coefficient) in KC.kinematic_factor(sector)
                atom_sector = KC._kinematic_atom_sector(sector, kinematic_monomial)
                base_sector = untransformed_sector(atom_sector)
                reflected_sector, reflected_support_factor =
                    KC._transform_kernel_sector(atom_sector, reflection)
                reflected_monomial = KC.transform_loop_momenta(
                    occupation_monomial, reflection
                )

                choose_reflected = candidate_isless(
                    reflected_sector,
                    reflected_monomial,
                    base_sector,
                    occupation_monomial,
                )
                selected_sector = choose_reflected ? reflected_sector : base_sector
                selected_monomial = choose_reflected ? reflected_monomial : occupation_monomial
                support_factor =
                    choose_reflected ? reflected_support_factor : one(Rational{Int})
                coefficient =
                    convert(D, occupation_coefficient) *
                    convert(D, kinematic_coefficient) *
                    convert(D, support_factor)
                KC._push_kernel_polynomial!(
                    out,
                    selected_sector,
                    KC.OccupationPolynomial{D,S}([selected_monomial => coefficient]),
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

function reduced_sector(sector::KC.CollisionKernelSector{S}) where {S<:KC.Statistics}
    return KC.ReducedCollisionSector{S}(
        KC.parameters(sector),
        KC.momentum_basis(sector),
        KC.external_wigner_momentum(sector),
        KC.kinematic_factor(sector),
        KC.frequency_support(sector),
    )
end

function gc_workspace(sector::KC.ReducedCollisionSector, monomial::KC.OccupationMonomial)
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

function gc_requotient_terms(
    expression::KC.LoopQuotientedExpression{C,S}
) where {C<:Number,S<:KC.Statistics}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{D,S}}()
    for (sector, polynomial) in KC.loop_quotient_terms(expression)
        source_sector = reduced_sector(sector)
        for (monomial, coefficient) in polynomial
            workspace = gc_workspace(source_sector, monomial)
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

function cases()
    out = Pair{String,Any}[]
    for (S, family) in ((Boson, direct_loop_ϕ), (Fermion, direct_loop_ψ))
        statistics = S === Boson ? "boson" : "fermion"
        for support in SUPPORT_KINDS,
            kinematic in KINEMATIC_KINDS,
            occupation in OCCUPATION_KINDS
            label = "$statistics/$support/$kinematic/$occupation"
            push!(
                out,
                label => one_loop_expression(S, family, support, kinematic, occupation),
            )
        end
    end
    return out
end

all_cases = cases()
semantic_failures = String[]
time_regressions = Tuple{String,Float64}[]
memory_regressions = Tuple{String,Int,Int}[]
worst_ratio = ("", 0.0)

@testset "direct one-loop semantics" begin
    for (label, expression) in all_cases
        direct = direct_one_loop_quotient(expression)
        nauty = nauty_quotient(expression)
        gc = KC.quotient_loop_momenta(expression)
        oracle = KC.loop_quotient_terms(gc)
        exact = gc_requotient_terms(direct) == oracle && gc_requotient_terms(nauty) == oracle
        exact || push!(semantic_failures, label)
        @test exact
    end
end

@testset "direct one-loop performance" begin
    global worst_ratio
    for (label, expression) in all_cases
        direct_one_loop_quotient(expression)
        nauty_quotient(expression)
        direct_trial = @benchmark direct_one_loop_quotient($expression) samples = 31 evals = 1
        nauty_trial = @benchmark nauty_quotient($expression) samples = 31 evals = 1
        direct = median(direct_trial)
        nauty = median(nauty_trial)
        if direct.time > nauty.time || direct.memory > nauty.memory
            direct_trial = @benchmark direct_one_loop_quotient($expression) samples = 81 evals = 1
            nauty_trial = @benchmark nauty_quotient($expression) samples = 81 evals = 1
            direct = median(direct_trial)
            nauty = median(nauty_trial)
        end
        ratio = direct.time / nauty.time
        ratio > worst_ratio[2] && (worst_ratio = (label, ratio))
        direct.time > nauty.time && push!(time_regressions, (label, ratio))
        direct.memory > nauty.memory &&
            push!(memory_regressions, (label, direct.memory, nauty.memory))
    end
    println("semantic failures: ", length(semantic_failures))
    println("direct one-loop time regressions: ", length(time_regressions))
    println("direct one-loop memory regressions: ", length(memory_regressions))
    println(
        "worst direct/Nauty ratio: ",
        round(worst_ratio[2]; digits=3),
        "x at ",
        worst_ratio[1],
    )
    @test isempty(time_regressions)
    @test isempty(memory_regressions)
end
