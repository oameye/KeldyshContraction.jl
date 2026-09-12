import KeldyshContraction as KC

@qfields benchmark_pwave_ψ::Fermion

function benchmark_fermionic_pwave_lagrangian()
    ψ₁ = benchmark_pwave_ψ[One]
    ψ₂ = benchmark_pwave_ψ[Two]
    bψ₁ = bar(ψ₁)
    bψ₂ = bar(ψ₂)
    ∂bψ₁ = partial(bψ₁, :x)
    ∂bψ₂ = partial(bψ₂, :x)
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus

    ψplus_minus = ψ₁(minus) + ψ₂(minus)
    ∂ψplus_minus = partial(ψ₁(minus), :x) + partial(ψ₂(minus), :x)
    ψplus_plus = ψ₁(plus) + ψ₂(plus)
    ∂ψplus_plus = partial(ψ₁(plus), :x) + partial(ψ₂(plus), :x)
    ψminus_plus = ψ₁(plus) - ψ₂(plus)
    ∂ψminus_plus = partial(ψ₁(plus), :x) - partial(ψ₂(plus), :x)

    bψplus = bψ₁ + bψ₂
    bψminus = bψ₂ - bψ₁
    ∂bψplus = ∂bψ₁ + ∂bψ₂
    ∂bψminus = ∂bψ₂ - ∂bψ₁

    Jplus_minus = ψplus_minus * ∂ψplus_minus
    Jplus_plus = ψplus_plus * ∂ψplus_plus
    Jminus_plus = ψminus_plus * ∂ψminus_plus
    Jplus_dagger = ∂bψplus * bψplus
    Jminus_dagger = ∂bψminus * bψminus

    loss =
        (1 // 8) *
        im *
        (
            (Jplus_dagger - Jminus_dagger) * Jplus_minus -
            (Jplus_plus - Jminus_plus) * Jminus_dagger
        )
    return InteractionLagrangian(loss, :γp)
end

function benchmark_fermionic_pwave_fixtures()
    L = benchmark_fermionic_pwave_lagrangian()
    G = DressedPropagator(L, Val(1), Val(3); simplify=true, _set_reg_to_zero=false)
    GF = fourier_transform(G)
    ΣF = SelfEnergy(GF)
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    kinetic = kinetic_expression(ΣW)
    off_shell = off_shell_collision_expression(kinetic)
    collision = spectral_dispersive_collision(off_shell)
    reduced = reduce_frequency_collision(collision)
    occupation = occupation_reduced_expression(reduced)
    quotient = quotient_loop_momenta(occupation)
    return L, G, GF, ΣF, ΣW, kinetic, off_shell, collision, reduced, occupation, quotient
end

function benchmark_fermionic_pwave_loss!(suite)
    L, G, GF, ΣF, ΣW, kinetic, off_shell, collision, reduced, occupation, quotient = benchmark_fermionic_pwave_fixtures()
    group = suite["Fermionic p-wave loss"]

    group["InteractionLagrangian"] = @benchmarkable benchmark_fermionic_pwave_lagrangian() seconds =
        10 evals = 1
    group["DressedPropagator"] = @benchmarkable DressedPropagator(
        $L, Val(1), Val(3); simplify=true, _set_reg_to_zero=false
    ) seconds = 10 evals = 1
    group["Fourier"] = @benchmarkable fourier_transform($G) seconds = 10 evals = 1
    group["SelfEnergy"] = @benchmarkable SelfEnergy($GF) seconds = 10 evals = 1
    group["Wigner"] = @benchmarkable wigner_transform($ΣF; gradient_order=Val(0)) seconds =
        10 evals = 1
    group["KineticSelfEnergy"] = @benchmarkable kinetic_expression($ΣW) seconds = 10 evals =
        1
    group["exact KB collision"] = @benchmarkable off_shell_collision_expression($kinetic) seconds =
        10 evals = 1
    group["spectral-dispersive"] = @benchmarkable spectral_dispersive_collision($off_shell) seconds =
        10 evals = 1
    group["frequency reduction"] = @benchmarkable reduce_frequency_collision($collision) seconds =
        10 evals = 1
    group["F to n"] = @benchmarkable occupation_reduced_expression($reduced) seconds = 10 evals =
        1
    group["loop quotient"] = @benchmarkable quotient_loop_momenta($occupation) seconds = 10 evals =
        1
    group["CollisionKernel"] = @benchmarkable collision_kernel($quotient) seconds = 10 evals =
        1
    group["complete pipeline"] = @benchmarkable begin
        Gp = DressedPropagator($L, Val(1), Val(3); simplify=true, _set_reg_to_zero=false)
        GFp = fourier_transform(Gp)
        ΣFp = SelfEnergy(GFp)
        ΣWp = wigner_transform(ΣFp; gradient_order=Val(0))
        kineticp = kinetic_expression(ΣWp)
        off_shellp = off_shell_collision_expression(kineticp)
        collisionp = spectral_dispersive_collision(off_shellp)
        reducedp = reduce_frequency_collision(collisionp)
        occupationp = occupation_reduced_expression(reducedp)
        quotientp = quotient_loop_momenta(occupationp)
        collision_kernel(quotientp)
    end seconds = 30 evals = 1
    return suite
end
