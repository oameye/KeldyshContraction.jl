using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields collision_api_ϕ::Boson
@qfields collision_api_ψ::Fermion

function collision_api_loss_propagator()
    c = collision_api_ϕ[Classical]
    q = collision_api_ϕ[Quantum]
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
    return DressedPropagator(L, Val(1), Val(3); preserve_regularisation=true)
end

function collision_api_pwave_loss_propagator()
    ψ1, ψ2 = collision_api_ψ[One], collision_api_ψ[Two]
    bψ1, bψ2 = bar(ψ1), bar(ψ2)

    ψplus = ψ1 + ψ2
    ∂ψplus = partial(ψ1, :x) + partial(ψ2, :x)
    ψminus = ψ1 - ψ2
    ∂ψminus = partial(ψ1, :x) - partial(ψ2, :x)
    bψplus = bψ1 + bψ2
    ∂bψplus = partial(bψ1, :x) + partial(bψ2, :x)
    bψminus = bψ2 - bψ1
    ∂bψminus = partial(bψ2, :x) - partial(bψ1, :x)

    Pplus_dagger = ∂bψplus * bψplus
    Pminus_dagger = ∂bψminus * bψminus

    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus
    ψplus_minus = ψ1(minus) + ψ2(minus)
    ∂ψplus_minus = partial(ψ1(minus), :x) + partial(ψ2(minus), :x)
    ψplus_plus = ψ1(plus) + ψ2(plus)
    ∂ψplus_plus = partial(ψ1(plus), :x) + partial(ψ2(plus), :x)
    ψminus_plus = ψ1(plus) - ψ2(plus)
    ∂ψminus_plus = partial(ψ1(plus), :x) - partial(ψ2(plus), :x)

    Pplus_minus = ψplus_minus * ∂ψplus_minus
    Pplus_plus = ψplus_plus * ∂ψplus_plus
    Pminus_plus = ψminus_plus * ∂ψminus_plus

    loss =
        (1 // 8) *
        im *
        (
            (Pplus_dagger - Pminus_dagger) * Pplus_minus -
            (Pplus_plus - Pminus_plus) * Pminus_dagger
        )
    L = InteractionLagrangian(loss, :γp)
    return DressedPropagator(L, Val(1), Val(3); preserve_regularisation=true)
end

function explicit_collision_pipeline(G)
    GF = fourier_transform(G)
    ΣF = SelfEnergy(GF)
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    KΣ = kinetic_expression(ΣW)
    I = off_shell_collision_expression(KΣ)
    SD = spectral_dispersive_collision(I)
    R = reduce_frequency_collision(SD)
    N = occupation_reduced_expression(R)
    Q = quotient_loop_momenta(N)
    C = collision_kernel(Q)
    return (; GF, ΣF, ΣW, KΣ, I, SD, R, N, Q, C)
end

@testset "high-level collision kernel compiler" begin
    G = collision_api_loss_propagator()
    stages = explicit_collision_pipeline(G)
    automatic = @inferred collision_kernel(G)
    reduced_direct = @inferred collision_kernel(stages.R)

    @test collision_kernel_terms(automatic) == collision_kernel_terms(stages.C)
    @test collision_kernel_terms(reduced_direct) == collision_kernel_terms(stages.C)
    @test target_family(automatic) === target_family(stages.C)
    @test parameters(automatic) == parameters(stages.C)
    @test KC.order(automatic) == KC.order(stages.C)
    @test KC.statistics(automatic) === KC.statistics(stages.C)
    @test isempty(reduced_blocked_terms(stages.R))
    @test isempty(reduced_trotter_terms(stages.R))
end

@testset "kinetic compiler displays are physics-facing" begin
    G = collision_api_loss_propagator()
    stages = explicit_collision_pipeline(G)

    displays = (
        (G, "G^K\\!\\left(x_1,x_2\\right)"),
        (stages.GF, "G_F^K\\!\\left(k\\right)"),
        (stages.ΣF, "\\Sigma_F^K\\!\\left(k\\right)"),
        (stages.ΣW, "\\Sigma_W^K\\!\\left(X,k\\right)"),
        (stages.KΣ, "\\Sigma_{\\mathrm{kin}}^K\\!\\left(k\\right)"),
        (stages.I, "I_{\\mathrm{coll}}(k)"),
        (stages.SD, "I_{\\mathrm{coll}}(k)"),
        (stages.R, "I_{\\mathrm{reg}}(k)"),
        (stages.N, "C_n^{\\mathrm{reg}}(k)"),
        (stages.Q, "C_n^{\\mathrm{quot}}(k)"),
        (stages.C, "C_n(k)"),
    )

    for (value, expected) in displays
        rendered = sprint(show, MIME"text/latex"(), value)
        @test startswith(rendered, "\$\$")
        @test occursin(expected, rendered)
        @test occursin("\\gamma", rendered)
        @test endswith(rendered, "\$\$")
        @test !occursin("\\[", rendered)
        @test !occursin(r"G\^[KRA]_\{[^}]+\}\^\{\[\\Delta t=", rendered)
    end

    routed = sprint(show, MIME"text/latex"(), stages.GF)
    @test occursin(r"G\^\{[KRA],\[\\Delta t=[+-]1\\,0\^\+\]\}_\{", routed)

    @test occursin("Kadanoff-Baym", sprint(show, MIME"text/plain"(), stages.I))
    @test occursin("causal blockers: 0", sprint(show, MIME"text/plain"(), stages.R))
    @test occursin("C_n(k) =", sprint(show, MIME"text/plain"(), stages.C))
end

@testset "collision-stage displays render computed p-wave physics" begin
    stages = explicit_collision_pipeline(collision_api_pwave_loss_propagator())

    sd = sprint(show, MIME"text/latex"(), stages.SD)
    reduced = sprint(show, MIME"text/latex"(), stages.R)
    occupation = sprint(show, MIME"text/latex"(), stages.N)
    quotiented = sprint(show, MIME"text/latex"(), stages.Q)
    kernel = sprint(show, MIME"text/latex"(), stages.C)

    @test occursin("γp", sd)
    @test occursin("A_{", sd)
    @test occursin("F_{", sd)
    @test occursin("_{x}", sd)

    @test occursin("γp", reduced)
    @test occursin("F_{", reduced)
    @test occursin("_{x}", reduced)

    for rendered in (occupation, quotiented, kernel)
        @test occursin("γp", rendered)
        @test occursin("n_{", rendered)
        @test occursin("_{x}", rendered)
        @test occursin("k", rendered)
        @test occursin("q", rendered)
    end

    @test !occursin("q_i\\sim", quotiented)
    @test !occursin("\\mathcal K", kernel)

    plain_kernel = sprint(show, MIME"text/plain"(), stages.C)
    @test occursin("γp", plain_kernel)
    @test occursin("n_", plain_kernel)
    @test occursin("_x", plain_kernel)
end