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
        :MomentumBasis,
        :matrix,
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
