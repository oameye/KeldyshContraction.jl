using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields quotient_pv_ϕ::Boson

@testset "principal-value orientation factor survives loop quotient" begin
    C = KC.ComplexRationals
    basis = KC.MomentumBasis(3)
    external = basis[1]
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    parameter = KC.ParameterMonomial(:g)

    energy(momentum) = KC.EnergyForm(KC.DispersionAtom(quotient_pv_ϕ, momentum))
    pv, pv_factor = KC.principal_value_support(energy(r) - energy(q))
    @test pv_factor == 1
    support = KC.FrequencySupport(KC.EnergyShell{Boson}[], [pv])
    kinematic = KC.MomentumPolynomial(KC.MomentumMonomial(), one(C))
    sector = KC.ReducedCollisionSector{Boson}(
        parameter, basis, external, kinematic, support
    )

    nq = OccupationPolynomial(OccupationAtom(quotient_pv_ϕ, q), one(C))
    expression = KC.OccupationReducedExpression{C,Boson,1,0,KC.HomogeneousWignerContext}(
        Dict(sector => nq), quotient_pv_ϕ, parameter, KC.HomogeneousWignerContext()
    )
    kernel = @inferred collision_kernel(expression)

    swap = loop_permutation_transform(basis, external, [2, 1])
    swapped_support, support_factor = KC._transform_frequency_support(support, swap)
    @test support_factor == -1
    swapped_sector = KC.ReducedCollisionSector{Boson}(
        parameter, basis, external, transform_loop_momenta(kinematic, swap), swapped_support
    )
    swapped_polynomial = support_factor * transform_loop_momenta(nq, swap)
    swapped_expression = KC.OccupationReducedExpression{
        C,Boson,1,0,KC.HomogeneousWignerContext
    }(
        Dict(swapped_sector => swapped_polynomial),
        quotient_pv_ϕ,
        parameter,
        KC.HomogeneousWignerContext(),
    )

    swapped_kernel = @inferred collision_kernel(swapped_expression)
    @test collision_kernel_terms(swapped_kernel) == collision_kernel_terms(kernel)
    @test all(
        key ->
            isempty(frequency_support(key).shells) &&
            length(frequency_support(key).principal_values) == 1,
        keys(collision_kernel_terms(kernel)),
    )
end
