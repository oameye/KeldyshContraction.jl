using KeldyshContraction, Test
using Combinatorics: Combinatorics
import KeldyshContraction as KC

@qfields projective_loop_ϕ::Boson

function projective_four_loop_fixture()
    C = KC.ComplexRationals
    basis = KC.MomentumBasis(5)
    external = basis[1]
    k = KC.basis_momentum(basis, 1)
    loops = [KC.basis_momentum(basis, index) for index in 2:5]
    q, r = loops[1], loops[2]
    p = -k + q + r

    energy(momentum) = KC.EnergyForm(KC.DispersionAtom(projective_loop_ϕ, momentum))
    shell, factor = KC.energy_shell(energy(k) + energy(p) - energy(q) - energy(r))
    @test factor == 1
    support = KC.FrequencySupport([shell], KC.PrincipalValueSupport{Boson}[])

    qx = KC.MomentumComponent(q, :x)
    rx = KC.MomentumComponent(r, :x)
    kinematic = KC.MomentumPolynomial(KC.MomentumMonomial([qx, rx]), one(C))
    parameter = KC.ParameterMonomial(:g)^2
    sector = KC.ReducedCollisionSector{Boson}(
        parameter, basis, external, kinematic, support
    )

    n(momentum) = OccupationPolynomial(OccupationAtom(projective_loop_ϕ, momentum), one(C))
    nk, np, nq, nr = n(k), n(p), n(q), n(r)
    polynomial = -4 * nk * np * nq - 2 * nk * np + 2 * nk * nq * nr + 2 * nq * nr + nr - nq
    for loop in Iterators.drop(loops, 2)
        polynomial += n(loop) + n(q + loop)
    end

    return KC.OccupationReducedExpression{C,Boson,2,0,KC.HomogeneousWignerContext}(
        Dict(sector => polynomial),
        projective_loop_ϕ,
        parameter,
        KC.HomogeneousWignerContext(),
    )
end

function signed_permutation_expression(
    expression::KC.OccupationReducedExpression{C,S,O,G,Ctx},
    permutation::AbstractVector{<:Integer},
    signs::AbstractVector{<:Integer},
) where {C<:Number,S<:KC.Statistics,O,G,Ctx<:KC.AbstractWignerContext}
    sector, polynomial = only(KC.occupation_reduced_terms(expression))
    basis = KC.momentum_basis(sector)
    external = KC.external_wigner_momentum(sector)
    transform = KC.loop_permutation_transform(
        basis, external, permutation, signs, zeros(Int, length(permutation))
    )
    support, support_factor = KC._transform_frequency_support(
        KC.frequency_support(sector), transform
    )
    transformed_sector = KC.ReducedCollisionSector{S}(
        KC.parameters(sector),
        basis,
        external,
        KC.transform_loop_momenta(KC.kinematic_factor(sector), transform),
        support,
    )
    transformed_polynomial =
        support_factor * KC.transform_loop_momenta(polynomial, transform)
    return KC.OccupationReducedExpression{C,S,O,G,Ctx}(
        Dict(transformed_sector => transformed_polynomial),
        KC.target_family(expression),
        KC.parameters(expression),
        KC.wigner_context(expression),
    )
end

@testset "projective derivative quotient is invariant under all signed loop permutations" begin
    expression = projective_four_loop_fixture()
    reference = KC.loop_quotient_terms(KC.quotient_loop_momenta(expression))
    checked = 0

    for permutation in Combinatorics.permutations(collect(1:4))
        for mask in 0:15
            signs = Int[isodd(mask >> (slot - 1)) ? -1 : 1 for slot in 1:4]
            transformed = signed_permutation_expression(expression, permutation, signs)
            quotient = KC.loop_quotient_terms(KC.quotient_loop_momenta(transformed))
            @test quotient == reference
            checked += 1
        end
    end

    @test checked == factorial(4) * 2^4
end
