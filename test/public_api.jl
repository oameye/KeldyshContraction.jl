using KeldyshContraction, Test
import KeldyshContraction as KC

function is_exported(name)
    return if VERSION >= v"1.11"
        Base.isexported(KC, name)
    else
        name in names(KC; all=false, imported=false)
    end
end

@testset "declared public API" begin
    # Ordinary workflow names are imported by `using KeldyshContraction`.
    for name in (
        :InteractionLagrangian,
        :DressedPropagator,
        :SelfEnergy,
        :fourier_transform,
        :wigner_transform,
        :kinetic_expression,
        :off_shell_collision_expression,
        :spectral_dispersive_collision,
        :reduce_frequency_collision,
        :occupation_reduced_expression,
        :quotient_loop_momenta,
        :collision_kernel,
    )
        @test is_exported(name)
    end

    # Qualified public names are semantic inspection used by the manual and canonical examples.
    for name in (
        :Regularisation,
        :field_families,
        :target_family,
        :parameters,
        :matrix,
        :order,
        :statistics,
        :keldysh_component,
        :retarded_component,
        :advanced_component,
        :kinematic_factor,
        :frequency_support,
        :reduced_regular_terms,
        :reduced_blocked_terms,
        :reduced_trotter_terms,
        :collision_kernel_terms,
    )
        @test !is_exported(name)
        VERSION >= v"1.11" && @test Base.ispublic(KC, name)
    end

    # Returned representation types and compiler machinery are implementation details.
    for name in (
        :QTerm,
        :QMul,
        :QAdd,
        :Edge,
        :Diagram,
        :Diagrams,
        :DressedPropagatorSum,
        :SelfEnergySum,
        :MomentumBasis,
        :FourierDressedPropagator,
        :FourierSelfEnergy,
        :WignerDressedPropagator,
        :WignerSelfEnergy,
        :KineticSelfEnergy,
        :OffShellCollisionExpression,
        :SpectralDispersiveCollision,
        :FrequencySupport,
        :ReducedFrequencyCollision,
        :OccupationReducedExpression,
        :LoopQuotientedExpression,
        :CollisionKernelSector,
        :KineticLine,
        :KineticTerm,
        :EnergyForm,
        :FullRankFrequencyReduction,
        :ExceptionalFrequencyClassification,
        :StatisticalPolynomial,
        :OccupationPolynomial,
        :LoopMomentumTransform,
        :CanonicalFrequencyCollision,
        :CausalFrequencyExpression,
        :TrotterFrequencyState,
    )
        @test !is_exported(name)
        VERSION >= v"1.11" && @test !Base.ispublic(KC, name)
    end

    # The repaired reducer has one typed blocker branch, not the historical split API.
    @test !isdefined(KC, :reduced_dependent_terms)
    @test !isdefined(KC, :reduced_causal_terms)
    @test !isdefined(KC, :Destroy)
    @test !isdefined(KC, :Create)
end

function test_rak_accessors(result)
    @test KC.keldysh_component(result) === result.keldysh
    @test KC.retarded_component(result) === result.retarded
    @test KC.advanced_component(result) === result.advanced
    return nothing
end

@testset "stable RAK accessors" begin
    @qfields api_ϕ::Boson
    c, q = api_ϕ[Classical], api_ϕ[Quantum]
    interaction = -(
        (1 // 2) * (c^2 + q^2) * bar(c) * bar(q) + (1 // 2) * c * q * (bar(c)^2 + bar(q)^2)
    )
    L = InteractionLagrangian(interaction, :g)
    G = DressedPropagator(L, Val(1), Val(3))
    Σ = SelfEnergy(G)
    GF = fourier_transform(G)
    ΣF = SelfEnergy(GF)
    GW = wigner_transform(GF; gradient_order=Val(0))
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    kinetic = kinetic_expression(ΣW)

    test_rak_accessors(G)
    test_rak_accessors(Σ)
    test_rak_accessors(GF)
    test_rak_accessors(ΣF)
    test_rak_accessors(GW)
    test_rak_accessors(ΣW)
    test_rak_accessors(kinetic)
end
