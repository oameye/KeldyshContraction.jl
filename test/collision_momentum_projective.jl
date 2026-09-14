using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields projective_quotient_ϕ::Boson

@testset "projective kinematic scalar survives loop quotient" begin
    C = KC.ComplexRationals
    basis = KC.MomentumBasis(2)
    external = basis[1]
    q = KC.basis_momentum(basis, 2)
    parameter = KC.ParameterMonomial(:λ)
    support = KC.FrequencySupport(
        KC.EnergyShell{Boson}[], KC.PrincipalValueSupport{Boson}[]
    )

    qx = KC.MomentumComponent(q, :x)
    minus_qx = KC.MomentumComponent(-q, :x)
    plus_kinematic = KC.MomentumPolynomial(KC.MomentumMonomial([qx]), one(C))
    minus_kinematic = KC.MomentumPolynomial(KC.MomentumMonomial([minus_qx]), one(C))

    minus_monomial, minus_coefficient = only(minus_kinematic)
    @test minus_monomial == KC.MomentumMonomial([qx])
    @test minus_coefficient == -1

    plus_sector = KC.ReducedCollisionSector{Boson}(
        parameter, basis, external, plus_kinematic, support
    )
    minus_sector = KC.ReducedCollisionSector{Boson}(
        parameter, basis, external, minus_kinematic, support
    )

    nq = OccupationPolynomial(OccupationAtom(projective_quotient_ϕ, q), one(C))
    expression = KC.OccupationReducedExpression{C,Boson,1,0,KC.HomogeneousWignerContext}(
        Dict(plus_sector => nq, minus_sector => nq),
        projective_quotient_ϕ,
        parameter,
        KC.HomogeneousWignerContext(),
    )

    sign_flip = loop_permutation_transform(basis, external, [1], [-1], [0])
    transformed = transform_loop_momenta(plus_kinematic, sign_flip)
    transformed_monomial, transformed_coefficient = only(transformed)
    @test transformed_monomial == KC.MomentumMonomial([qx])
    @test transformed_coefficient == -1

    kernel = @inferred collision_kernel(expression)
    @test isempty(collision_kernel_terms(kernel))
end
