using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields quotient_test_ϕ::Boson
@qfields quotient_test_ψ::Fermion

function quotient_identity_sector(
    ::Type{S}, family, basis, external, parameter
) where {S<:KC.Statistics}
    C = KC.ComplexRationals
    kinematic = KC.MomentumPolynomial(KC.MomentumMonomial(), one(C))
    support = KC.FrequencySupport(KC.EnergyShell{S}[], KC.PrincipalValueSupport{S}[])
    return KC.ReducedCollisionSector{S}(parameter, basis, external, kinematic, support)
end

function quotient_transform_sector(
    sector::KC.ReducedCollisionSector{S}, transform::LoopMomentumTransform
) where {S<:KC.Statistics}
    support, factor = KC._transform_frequency_support(frequency_support(sector), transform)
    transformed = KC.ReducedCollisionSector{S}(
        parameters(sector),
        KC.momentum_basis(sector),
        external_wigner_momentum(sector),
        transform_loop_momenta(KC.kinematic_factor(sector), transform),
        support,
    )
    return transformed, factor
end

function quotient_expression(
    sector::KC.ReducedCollisionSector{S},
    polynomial::OccupationPolynomial{C,S},
    family::FieldFamily{S},
    parameter::KC.ParameterMonomial,
    ::Val{O}=Val(2),
) where {C<:Number,S<:KC.Statistics,O}
    return KC.OccupationReducedExpression{C,S,O,0,KC.HomogeneousWignerContext}(
        Dict(sector => polynomial), family, parameter, KC.HomogeneousWignerContext()
    )
end

@testset "exact loop-momentum transforms" begin
    basis = KC.MomentumBasis(3)
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    external = basis[1]

    identity = loop_permutation_transform(basis, external, [1, 2])
    @test transform_loop_momenta(k, identity) == k
    @test transform_loop_momenta(q, identity) == q
    @test transform_loop_momenta(r, identity) == r

    swap = loop_permutation_transform(basis, external, [2, 1])
    @test transform_loop_momenta(k, swap) == k
    @test transform_loop_momenta(q, swap) == r
    @test transform_loop_momenta(r, swap) == q
    @test transform_loop_momenta(-k + q + r, swap) == -k + q + r
    @test loop_transform_matrix(swap) == Rational{Int}[1 0 0; 0 0 1; 0 1 0]

    sign_flip = loop_permutation_transform(basis, external, [1, 2], [-1, 1], [0, 0])
    @test transform_loop_momenta(k, sign_flip) == k
    @test transform_loop_momenta(q, sign_flip) == -q
    @test transform_loop_momenta(r, sign_flip) == r
    @test loop_transform_matrix(sign_flip) == Rational{Int}[1 0 0; 0 -1 0; 0 0 1]

    reflected_shifted = loop_permutation_transform(basis, external, [1, 2], [-1, 1], [1, 0])
    @test transform_loop_momenta(q, reflected_shifted) == k - q
    @test transform_loop_momenta(r, reflected_shifted) == r

    nonunimodular = Rational{Int}[1 0 0; 0 2 0; 0 0 1]
    @test_throws ArgumentError LoopMomentumTransform(nonunimodular, 1)
    @test_throws ArgumentError LoopMomentumTransform(nonunimodular, Int16(1))

    moves_external = Rational{Int}[1 1 0; 0 1 0; 0 0 1]
    @test_throws ArgumentError LoopMomentumTransform(moves_external, 1)

    owned = Rational{Int}[1 0 0; 0 1 0; 0 0 1]
    transform = LoopMomentumTransform(owned, 1)
    owned[2, 2] = 2
    @test loop_transform_matrix(transform) == Rational{Int}[1 0 0; 0 1 0; 0 0 1]
end

@testset "termwise signed-permutation quotient" begin
    C = KC.ComplexRationals
    basis = KC.MomentumBasis(3)
    external = basis[1]
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    parameter = KC.ParameterMonomial(:g)
    sector = quotient_identity_sector(Boson, quotient_test_ϕ, basis, external, parameter)

    n(momentum) = OccupationPolynomial(OccupationAtom(quotient_test_ϕ, momentum), one(C))
    raw = n(q) + 3 * n(r)
    expression = quotient_expression(sector, raw, quotient_test_ϕ, parameter)

    quotient = @inferred quotient_loop_momenta(expression)
    kernel = @inferred collision_kernel(quotient)
    @test collision_kernel(expression).terms == kernel.terms
    @test length(loop_quotient_terms(quotient)) == 1
    _, polynomial = only(loop_quotient_terms(quotient))
    @test length(polynomial) == 1
    _, coefficient = only(polynomial)
    @test coefficient == 4

    swap = loop_permutation_transform(basis, external, [2, 1])
    swapped = quotient_expression(
        sector, transform_loop_momenta(raw, swap), quotient_test_ϕ, parameter
    )
    @test collision_kernel(swapped).terms == kernel.terms

    sign_flip = loop_permutation_transform(basis, external, [1, 2], [-1, 1], [0, 0])
    flipped = quotient_expression(
        sector, transform_loop_momenta(raw, sign_flip), quotient_test_ϕ, parameter
    )
    @test collision_kernel(flipped).terms == kernel.terms

    inequivalent = n(q) + n(q + r)
    inequivalent_expression = quotient_expression(
        sector, inequivalent, quotient_test_ϕ, parameter
    )
    _, inequivalent_polynomial = only(
        collision_kernel_terms(collision_kernel(inequivalent_expression))
    )
    @test length(inequivalent_polynomial) == 2
end

@testset "support and derivative kinematics canonicalize coherently" begin
    C = KC.ComplexRationals
    basis = KC.MomentumBasis(3)
    external = basis[1]
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    parameter = KC.ParameterMonomial(:λ)

    energy(momentum) = KC.EnergyForm(KC.DispersionAtom(quotient_test_ϕ, momentum))
    shell, shell_factor = KC.energy_shell(energy(k) + energy(q) - energy(r))
    @test shell_factor == 1
    support = KC.FrequencySupport([shell], KC.PrincipalValueSupport{Boson}[])
    qx = KC.MomentumComponent(q, :x)
    kinematic = KC.MomentumPolynomial(KC.MomentumMonomial([qx]), one(C))
    sector = KC.ReducedCollisionSector{Boson}(
        parameter, basis, external, kinematic, support
    )
    nq = OccupationPolynomial(OccupationAtom(quotient_test_ϕ, q), one(C))
    expression = quotient_expression(sector, nq, quotient_test_ϕ, parameter, Val(1))
    kernel = @inferred collision_kernel(expression)

    swap = loop_permutation_transform(basis, external, [2, 1])
    swapped_sector, swapped_factor = quotient_transform_sector(sector, swap)
    @test swapped_factor == 1
    swapped_expression = quotient_expression(
        swapped_sector, transform_loop_momenta(nq, swap), quotient_test_ϕ, parameter, Val(1)
    )
    @test collision_kernel(swapped_expression).terms == kernel.terms

    sign_flip = loop_permutation_transform(basis, external, [1, 2], [-1, 1], [0, 0])
    flipped_sector, flipped_factor = quotient_transform_sector(sector, sign_flip)
    @test flipped_factor == 1
    flipped_expression = quotient_expression(
        flipped_sector,
        transform_loop_momenta(nq, sign_flip),
        quotient_test_ϕ,
        parameter,
        Val(1),
    )
    @test collision_kernel(flipped_expression).terms == kernel.terms
end

@testset "fermionic derivative kinematics use the same quotient" begin
    C = KC.ComplexRationals
    basis = KC.MomentumBasis(3)
    external = basis[1]
    q = KC.basis_momentum(basis, 2)
    r = KC.basis_momentum(basis, 3)
    parameter = KC.ParameterMonomial(:γp)

    qx = KC.MomentumComponent(q, :x)
    rx = KC.MomentumComponent(r, :x)
    p_wave = KC.MomentumPolynomial(KC.MomentumMonomial([qx, rx]), one(C))
    support = KC.FrequencySupport(
        KC.EnergyShell{Fermion}[], KC.PrincipalValueSupport{Fermion}[]
    )
    sector = KC.ReducedCollisionSector{Fermion}(parameter, basis, external, p_wave, support)
    n(momentum) = OccupationPolynomial(OccupationAtom(quotient_test_ψ, momentum), one(C))
    expression = quotient_expression(sector, n(q) + 3 * n(r), quotient_test_ψ, parameter)

    kernel = @inferred collision_kernel(expression)
    @test KC.statistics(kernel) === Fermion
    @test length(collision_kernel_terms(kernel)) == 1
    kernel_sector, polynomial = only(collision_kernel_terms(kernel))
    @test KC.statistics(kernel_sector) === Fermion
    @test !iszero(KC.kinematic_factor(kernel_sector))
    @test length(polynomial) == 1
    _, coefficient = only(polynomial)
    @test coefficient == 4

    swap = loop_permutation_transform(basis, external, [2, 1])
    swapped_sector, factor = quotient_transform_sector(sector, swap)
    @test factor == 1
    swapped_expression = quotient_expression(
        swapped_sector,
        transform_loop_momenta(n(q) + 3 * n(r), swap),
        quotient_test_ψ,
        parameter,
    )
    @test collision_kernel(swapped_expression).terms == kernel.terms
end

@testset "canonical quotient is not specialized to two loops" begin
    C = KC.ComplexRationals
    basis = KC.MomentumBasis(5)
    external = basis[1]
    loops = [KC.basis_momentum(basis, i) for i in 2:5]
    parameter = KC.ParameterMonomial(:u)
    sector = quotient_identity_sector(Boson, quotient_test_ϕ, basis, external, parameter)
    n(momentum) = OccupationPolynomial(OccupationAtom(quotient_test_ϕ, momentum), one(C))
    raw = n(loops[1] + loops[3]) + 2 * n(loops[2] - loops[4])
    expression = quotient_expression(sector, raw, quotient_test_ϕ, parameter)
    kernel = @inferred collision_kernel(expression)

    transform = loop_permutation_transform(
        basis, external, [4, 2, 1, 3], [-1, 1, -1, 1], [0, 0, 0, 0]
    )
    transformed_expression = quotient_expression(
        sector, transform_loop_momenta(raw, transform), quotient_test_ϕ, parameter
    )
    @test collision_kernel(transformed_expression).terms == kernel.terms
end
