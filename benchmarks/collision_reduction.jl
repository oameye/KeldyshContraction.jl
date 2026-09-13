import KeldyshContraction as KC

@qfields benchmark_collision_ϕ::Boson

function benchmark_collision_fixture()
    c = benchmark_collision_ϕ[Classical]
    q = benchmark_collision_ϕ[Quantum]
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus
    loss =
        (1 // 2) *
        im *
        (
            bar(c) * bar(q) * (c(minus)^2 + q(minus)^2) -
            c(plus) * q(plus) * (bar(c)^2 + bar(q)^2) +
            2 * bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
        )

    L = InteractionLagrangian(loss, :γ)
    G = DressedPropagator(L, Val(1), Val(3); simplify=true, _set_reg_to_zero=false)
    GF = fourier_transform(G)
    ΣW = wigner_transform(SelfEnergy(GF); gradient_order=Val(0))
    return spectral_dispersive_collision(
        off_shell_collision_expression(kinetic_expression(ΣW))
    )
end

function benchmark_statistical_terms()
    k = KC.LinearMomentum([1, 0])
    q = KC.LinearMomentum([0, 1])
    Fk = StatisticalAtom(benchmark_collision_ϕ, k)
    Fq = StatisticalAtom(benchmark_collision_ϕ, q)
    one_monomial = StatisticalMonomial{Boson}()
    k_monomial = StatisticalMonomial([Fk])
    q_monomial = StatisticalMonomial([Fq])
    kq_monomial = StatisticalMonomial([Fq, Fk])
    C = KC.ComplexRationals
    return Pair{StatisticalMonomial{Boson},C}[
        kq_monomial => C(-2),
        one_monomial => C(-2),
        q_monomial => C(2),
        k_monomial => C(2),
        kq_monomial => C(-2),
        kq_monomial => C(2),
    ]
end

function benchmark_loop_quotient_fixture(nloops::Int=2)
    C = KC.ComplexRationals
    basis = KC.MomentumBasis(nloops + 1)
    external = basis[1]
    k = KC.basis_momentum(basis, 1)
    loops = [KC.basis_momentum(basis, i) for i in 2:(nloops + 1)]
    q = loops[1]
    r = loops[min(2, nloops)]
    p = -k + q + r

    energy(momentum) = KC.EnergyForm(KC.DispersionAtom(benchmark_collision_ϕ, momentum))
    shell, factor = KC.energy_shell(energy(k) + energy(p) - energy(q) - energy(r))
    factor == 1 || error("benchmark shell unexpectedly changed normalization")
    support = KC.FrequencySupport([shell], KC.PrincipalValueSupport{Boson}[])

    qx = KC.MomentumComponent(q, :x)
    rx = KC.MomentumComponent(r, :x)
    kinematic = KC.MomentumPolynomial(KC.MomentumMonomial([qx, rx]), one(C))
    sector = KC.ReducedCollisionSector{Boson}(
        KC.ParameterMonomial(:g)^2, basis, external, kinematic, support
    )

    n(momentum) =
        OccupationPolynomial(OccupationAtom(benchmark_collision_ϕ, momentum), one(C))
    nk, np, nq, nr = n(k), n(p), n(q), n(r)
    raw = -4 * nk * np * nq - 2 * nk * np + 2 * nk * nq * nr + 2 * nq * nr + nr - nq
    for loop in Iterators.drop(loops, 2)
        raw += n(loop) + n(q + loop)
    end
    return KC.OccupationReducedExpression{C,Boson,2,0,KC.HomogeneousWignerContext}(
        Dict(sector => raw),
        benchmark_collision_ϕ,
        KC.ParameterMonomial(:g)^2,
        KC.HomogeneousWignerContext(),
    )
end

function benchmark_collision_reduction!(suite)
    collision = benchmark_collision_fixture()
    reduced = reduce_frequency_collision(collision)
    statistical = only(values(reduced_regular_terms(reduced)))
    raw_terms = benchmark_statistical_terms()
    occupation2 = benchmark_loop_quotient_fixture(2)
    occupation4 = benchmark_loop_quotient_fixture(4)
    quotient = quotient_loop_momenta(occupation2)

    suite["Collision reduction"]["collision-level assembly"] = @benchmarkable reduce_frequency_collision(
        $collision
    ) seconds = 10
    suite["Collision reduction"]["statistical canonicalization"] = @benchmarkable StatisticalPolynomial{
        KC.ComplexRationals,Boson
    }(
        copy($raw_terms)
    ) seconds = 10
    suite["Collision reduction"]["occupation expansion"] = @benchmarkable occupation_collision_polynomial(
        $statistical
    ) seconds = 10
    suite["Collision reduction"]["loop-momentum quotient 2-loop"] = @benchmarkable quotient_loop_momenta(
        $occupation2
    ) seconds = 10
    suite["Collision reduction"]["loop-momentum quotient 4-loop"] = @benchmarkable quotient_loop_momenta(
        $occupation4
    ) seconds = 10
    suite["Collision reduction"]["CollisionKernel lowering"] = @benchmarkable collision_kernel(
        $quotient
    ) seconds = 10
    return suite
end
