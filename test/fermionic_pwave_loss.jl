using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields pwave_trotter_ϕ::Boson pwave_loss_ψ::Fermion

function bosonic_branch_trotter_loss()
    c = pwave_trotter_ϕ[Classical]
    q = pwave_trotter_ϕ[Quantum]
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus

    Jplus_minus = (c(minus) + q(minus))^2
    Jplus_plus = (c(plus) + q(plus))^2
    Jminus_plus = (c(plus) - q(plus))^2
    Jplus_dagger = (bar(c) + bar(q))^2
    Jminus_dagger = (bar(c) - bar(q))^2
    return (1 // 8) *
           im *
           (
               (Jplus_dagger - Jminus_dagger) * Jplus_minus -
               (Jplus_plus - Jminus_plus) * Jminus_dagger
           )
end

function bosonic_known_trotter_loss()
    c = pwave_trotter_ϕ[Classical]
    q = pwave_trotter_ϕ[Quantum]
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus
    return (1 // 2) *
           im *
           (
               bar(c) * bar(q) * (c(minus)^2 + q(minus)^2) -
               c(plus) * q(plus) * (bar(c)^2 + bar(q)^2) +
               2 * bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
           )
end

function collected_quartic_coefficients(
    expression::KC.QAdd{C,S}
) where {C<:Number,S<:KC.Statistics}
    out = Dict{NTuple{4,KC.Field{S}},C}()
    for term in KC.terms(expression)
        fields = KC.fields(term)
        length(fields) == 4 || throw(ArgumentError("expected quartic field monomial"))
        key = ntuple(i -> fields[i], Val(4))
        combined = get(out, key, zero(C)) + KC.coefficient(term)
        if iszero(combined)
            delete!(out, key)
        else
            out[key] = combined
        end
    end
    return out
end

function fermionic_pwave_loss_lagrangian()
    ψ₁ = pwave_loss_ψ[One]
    ψ₂ = pwave_loss_ψ[Two]
    bψ₁ = bar(ψ₁)
    bψ₂ = bar(ψ₂)
    ∂bψ₁ = partial(bψ₁, :x)
    ∂bψ₂ = partial(bψ₂, :x)
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus

    # Asymmetric fermionic LO inverse rotation:
    #
    #   ψ+  = (ψ₁+ψ₂)/√2,      ψ-  = (ψ₁-ψ₂)/√2,
    #   ψ̄+ = (ψ̄₁+ψ̄₂)/√2,     ψ̄- = (ψ̄₂-ψ̄₁)/√2.
    #
    # The Trotter-ordered Lindblad contour form used by the bosonic loss vertex can be
    # written for an arbitrary even jump as
    #
    #   i/2 [(J+†-J-†) J+^(-) - (J+^(+)-J-^(+)) J-†].
    #
    # Removing the shifts gives
    #   -i [J-†J+ - 1/2 J+†J+ - 1/2 J-†J-].
    #
    # For J = ψ ∂xψ each branch jump contains two factors 1/√2, so using the
    # unnormalised LO sums below leaves the exact overall coefficient i/8.
    ψplus_minus = ψ₁(minus) + ψ₂(minus)
    ∂ψplus_minus = partial(ψ₁(minus), :x) + partial(ψ₂(minus), :x)
    ψplus_plus = ψ₁(plus) + ψ₂(plus)
    ∂ψplus_plus = partial(ψ₁(plus), :x) + partial(ψ₂(plus), :x)
    ψminus_plus = ψ₁(plus) - ψ₂(plus)
    ∂ψminus_plus = partial(ψ₁(plus), :x) - partial(ψ₂(plus), :x)

    bψplus = bψ₁ + bψ₂
    bψminus = bψ₂ - bψ₁
    ∂bψplus = ∂bψ₁ + ∂bψ₂
    ∂bψminus = ∂bψ₂ - ∂bψ₁

    Jplus_minus = ψplus_minus * ∂ψplus_minus
    Jplus_plus = ψplus_plus * ∂ψplus_plus
    Jminus_plus = ψminus_plus * ∂ψminus_plus
    Jplus_dagger = ∂bψplus * bψplus
    Jminus_dagger = ∂bψminus * bψminus

    loss =
        (1 // 8) *
        im *
        (
            (Jplus_dagger - Jminus_dagger) * Jplus_minus -
            (Jplus_plus - Jminus_plus) * Jminus_dagger
        )
    return InteractionLagrangian(loss, :γp)
end

function fermionic_pwave_generated_kernel()
    # Building the branch-expanded symbolic fixture is intentionally outside the solver
    # inference contract. From DressedPropagator onward every public computational stage
    # is required to infer concretely.
    L = fermionic_pwave_loss_lagrangian()
    G = @inferred DressedPropagator(
        L, Val(1), Val(3); simplify=true, _set_reg_to_zero=false
    )
    GF = @inferred fourier_transform(G)
    ΣF = @inferred SelfEnergy(GF)
    ΣW = @inferred wigner_transform(ΣF; gradient_order=Val(0))
    kinetic = @inferred kinetic_expression(ΣW)
    off_shell = @inferred off_shell_collision_expression(kinetic)
    collision = @inferred spectral_dispersive_collision(off_shell)
    reduced = @inferred reduce_frequency_collision(collision)
    occupation = @inferred occupation_reduced_expression(reduced)
    quotient = @inferred quotient_loop_momenta(occupation)
    kernel = @inferred collision_kernel(quotient)
    return L,
    G, GF, ΣF, ΣW, kinetic, off_shell, collision, reduced, occupation, quotient,
    kernel
end

function pwave_momentum_component(momentum)
    C = KC.ComplexRationals
    component = KC.MomentumComponent(momentum, :x)
    monomial = KC.MomentumMonomial(KC.MomentumComponent[component])
    return KC.MomentumPolynomial(monomial, one(C))
end

function expected_pwave_loss_kernel(occupation, sector)
    basis = KC.momentum_basis(sector)
    external = KC.external_wigner_momentum(sector)
    external_index = only(i for (i, variable) in enumerate(basis) if variable == external)
    loop_index = only(i for i in eachindex(basis.variables) if i != external_index)
    k = KC.basis_momentum(basis, external_index)
    q = KC.basis_momentum(basis, loop_index)

    relative = pwave_momentum_component(k) - pwave_momentum_component(q)
    relative_squared = relative * relative
    support = KC.FrequencySupport(
        KC.EnergyShell{Fermion}[], KC.PrincipalValueSupport{Fermion}[]
    )
    expected_sector = KC.ReducedCollisionSector{Fermion}(
        KC.ParameterMonomial(:γp), basis, external, relative_squared, support
    )

    C = KC.ComplexRationals
    n(momentum) = KC.OccupationPolynomial(KC.OccupationAtom(pwave_loss_ψ, momentum), one(C))
    expected_polynomial = -n(k) * n(q)
    expected = typeof(occupation)(
        Dict(expected_sector => expected_polynomial),
        KC.target_family(occupation),
        KC.parameters(occupation),
        KC.wigner_context(occupation),
    )
    return @inferred collision_kernel(expected)
end

@testset "branch-space Trotter calibration" begin
    @test collected_quartic_coefficients(bosonic_branch_trotter_loss()) ==
        collected_quartic_coefficients(bosonic_known_trotter_loss())
end

@testset "fermionic p-wave loss contour vertex" begin
    L = fermionic_pwave_loss_lagrangian()
    @test L isa InteractionLagrangian{KC.ComplexRationals,Fermion}
    @test KC.field_families(L) == [pwave_loss_ψ]
    @test KC.target_family(L) === pwave_loss_ψ
    @test KC.parameters(L) == KC.ParameterMonomial(:γp)
end

@testset "end-to-end fermionic p-wave two-body loss" begin
    L, G, GF, ΣF, ΣW, kinetic, off_shell, collision, reduced, occupation, quotient, kernel = fermionic_pwave_generated_kernel()

    shared = (
        G, GF, ΣF, ΣW, kinetic, off_shell, collision, reduced, occupation, quotient, kernel
    )
    @test all(x -> KC.statistics(x) === Fermion, shared)
    @test all(x -> KC.target_family(x) === pwave_loss_ψ, (L, shared...))
    @test all(x -> KC.parameters(x) == KC.ParameterMonomial(:γp), (L, shared...))
    @test KC.order(G) == 1
    @test KC.order(ΣF) == 1
    @test KC.order(ΣW) == 1
    @test KC.gradient_order(ΣW) == Val(0)
    @test KC.gradient_order(kernel) == Val(0)

    @test isempty(KC.reduced_dependent_terms(reduced))
    @test isempty(KC.reduced_causal_terms(reduced))
    @test isempty(KC.reduced_trotter_terms(reduced))
    regular = KC.occupation_reduced_terms(occupation)
    @test !isempty(regular)
    for (sector, _) in regular
        @test isempty(KC.frequency_support(sector).shells)
        @test isempty(KC.frequency_support(sector).principal_values)
    end

    sector = first(keys(regular))
    expected = expected_pwave_loss_kernel(occupation, sector)
    @test KC.collision_kernel_terms(kernel) == KC.collision_kernel_terms(expected)
    @test KC.parameters(kernel) == KC.ParameterMonomial(:γp)
end
