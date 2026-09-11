using KeldyshContraction, Test

const KC = KeldyshContraction

@testset "collision statistical algebra" begin
    ϕ = FieldFamily{Boson}(:ϕ)
    k = KC.LinearMomentum([1, 0])
    q = KC.LinearMomentum([0, 1])

    Fk_atom = StatisticalAtom(ϕ, k)
    Fq_atom = StatisticalAtom(ϕ, q)

    @test StatisticalMonomial([Fq_atom, Fk_atom]) == StatisticalMonomial([Fk_atom, Fq_atom])
    @test hash(StatisticalMonomial([Fq_atom, Fk_atom])) ==
        hash(StatisticalMonomial([Fk_atom, Fq_atom]))

    Fk = StatisticalPolynomial(Fk_atom, 1 // 1)
    Fq = StatisticalPolynomial(Fq_atom, 1 // 1)
    oneF = one(typeof(Fk))

    @test occupation_statistics_sign(Boson) == 1
    @test occupation_collision_factor(Boson) == 1 // 2
    @test statistical_occupation_coefficients(Boson) == (Int8(1), Int8(2))

    nk = OccupationPolynomial(OccupationAtom(ϕ, k), 1 // 1)
    nq = OccupationPolynomial(OccupationAtom(ϕ, q), 1 // 1)

    @test @inferred(occupation_substitute(Fk)) == one(typeof(nk)) + 2 * nk
    @test iszero(@inferred(occupation_collision_polynomial(zero(typeof(Fk)))))

    # Independent first-order two-body-loss oracle:
    # I_package/γ = 2(F_k - 1)(1 - F_q), C_n = I_package/2.
    loss_F = 2 * (Fk - oneF) * (oneF - Fq)
    @test @inferred(occupation_collision_polynomial(loss_F)) == -4 * nk * nq

    ψ = FieldFamily{Fermion}(:ψ)
    p = KC.LinearMomentum([1])
    Fp_atom = StatisticalAtom(ψ, p)
    Fp = StatisticalPolynomial(Fp_atom, 1 // 1)
    np = OccupationPolynomial(OccupationAtom(ψ, p), 1 // 1)

    @test occupation_statistics_sign(Fermion) == -1
    @test occupation_collision_factor(Fermion) == -1 // 2
    @test statistical_occupation_coefficients(Fermion) == (Int8(1), Int8(-2))
    @test @inferred(occupation_substitute(Fp)) == one(typeof(np)) - 2 * np

    # Synthetic fermion normalization oracle: 2(F_p-1) -> +2 n_p after C_n=-I/2.
    fermion_collision = 2 * (Fp - one(typeof(Fp)))
    @test @inferred(occupation_collision_polynomial(fermion_collision)) == 2 * np

    Fk_exact = StatisticalPolynomial(Fk_atom, KC.ComplexRationals(1))
    nk_exact = OccupationPolynomial(OccupationAtom(ϕ, k), KC.ComplexRationals(1))
    @test @inferred(occupation_substitute(Fk_exact)) ==
        one(typeof(nk_exact)) + KC.ComplexRationals(2) * nk_exact
    @test iszero(@inferred(occupation_collision_polynomial(zero(typeof(Fk_exact)))))

    Fk_float = StatisticalPolynomial(Fk_atom, ComplexF64(1))
    nk_float = OccupationPolynomial(OccupationAtom(ϕ, k), ComplexF64(1))
    @test @inferred(occupation_substitute(Fk_float)) ==
        one(typeof(nk_float)) + ComplexF64(2) * nk_float
    @test iszero(@inferred(occupation_collision_polynomial(zero(typeof(Fk_float)))))

    # The statistics map must not touch derivative/p-wave kinematics. Construct one synthetic
    # regular fermionic sector with an explicit routed momentum component and verify that the
    # exact sector identity survives occupation projection unchanged.
    fermion_basis = KC.MomentumBasis(1)
    routed_p = KC.basis_momentum(fermion_basis, 1)
    component = KC.MomentumComponent(routed_p, :x)
    p_wave = KC.MomentumPolynomial(
        KC.MomentumMonomial(KC.MomentumComponent[component]), one(KC.ComplexRationals)
    )
    empty_support = KC.FrequencySupport(
        KC.EnergyShell{Fermion}[], KC.PrincipalValueSupport{Fermion}[]
    )
    fermion_parameter = KC.ParameterMonomial(:λ)
    fermion_sector = ReducedCollisionSector{Fermion}(
        fermion_parameter, fermion_basis, fermion_basis[1], p_wave, empty_support
    )
    context = HomogeneousWignerContext()
    C = Rational{Int}
    reduced_fermion = ReducedFrequencyCollision{C,Fermion,0,0,0,0,HomogeneousWignerContext}(
        Dict(fermion_sector => fermion_collision),
        Dict{ReducedDependentCollisionSector{Fermion},StatisticalPolynomial{C,Fermion}}(),
        Dict{ReducedCausalCollisionSector{Fermion},StatisticalPolynomial{C,Fermion}}(),
        ReducedTrotterCollisionTerm{C,Fermion,0,0}[],
        ψ,
        fermion_parameter,
        context,
    )
    fermion_occupation = @inferred occupation_reduced_expression(reduced_fermion)
    occupation_sector, occupation_polynomial = only(
        occupation_reduced_terms(fermion_occupation)
    )
    @test occupation_sector == fermion_sector
    @test KC.kinematic_factor(occupation_sector) == p_wave
    @test occupation_polynomial == 2 * np
end
