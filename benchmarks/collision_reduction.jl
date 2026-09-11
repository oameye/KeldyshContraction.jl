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

function benchmark_collision_reduction!(suite)
    collision = benchmark_collision_fixture()
    reduced = reduce_frequency_collision(collision)
    statistical = only(values(reduced_regular_terms(reduced)))
    raw_terms = benchmark_statistical_terms()

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
    return suite
end
