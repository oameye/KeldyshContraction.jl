using KeldyshContraction, Test
import KeldyshContraction as KC

@testset "declared public API" begin
    exported = Set(names(KC; all=false, imported=false))

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
        @test name in exported
    end

    # Stable representation/inspection API is qualified, not namespace-polluting.
    for name in (
        :QTerm,
        :Regularisation,
        :Edge,
        :Diagram,
        :Diagrams,
        :contractions,
        :topology,
        :propagator_type,
        :DressedPropagatorSum,
        :matrix,
        :order,
        :statistics,
        :keldysh_component,
        :retarded_component,
        :advanced_component,
        :MomentumBasis,
        :FourierSelfEnergy,
        :WignerSelfEnergy,
        :KineticSelfEnergy,
        :OffShellCollisionExpression,
        :SpectralDispersiveCollision,
        :ReducedFrequencyCollision,
        :OccupationReducedExpression,
        :LoopQuotientedExpression,
        :CollisionKernelSector,
    )
        @test name ∉ exported
        VERSION >= v"1.11" && @test Base.ispublic(KC, name)
    end

    # Algorithmic IR and canonicalization machinery are implementation details.
    for name in (
        :KineticLine,
        :KineticTerm,
        :EnergyForm,
        :FullRankFrequencyReduction,
        :ExceptionalFrequencyClassification,
        :StatisticalPolynomial,
        :OccupationPolynomial,
        :LoopMomentumTransform,
    )
        @test name ∉ exported
        VERSION >= v"1.11" && @test !Base.ispublic(KC, name)
    end

    @test !isdefined(KC, :Destroy)
    @test !isdefined(KC, :Create)
end

@testset "stable RAK accessors" begin
    @qfields api_ϕ::Boson
    c, q = api_ϕ[Classical], api_ϕ[Quantum]
    interaction = -(
        (1 // 2) * (c^2 + q^2) * bar(c) * bar(q) +
        (1 // 2) * c * q * (bar(c)^2 + bar(q)^2)
    )
    L = InteractionLagrangian(interaction, :g)
    G = DressedPropagator(L, Val(1), Val(3))
    Σ = SelfEnergy(G)

    @test KC.keldysh_component(G) === G.keldysh
    @test KC.retarded_component(G) === G.retarded
    @test KC.advanced_component(G) === G.advanced
    @test KC.keldysh_component(Σ) === Σ.keldysh
    @test KC.retarded_component(Σ) === Σ.retarded
    @test KC.advanced_component(Σ) === Σ.advanced
end
