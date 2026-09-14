using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields collision_api_ϕ::Boson

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
        (G, "G^{R,A,K}"),
        (stages.GF, "G_F^{R,A,K}"),
        (stages.ΣF, "\\Sigma_F^{R,A,K}"),
        (stages.ΣW, "\\Sigma_W^{R,A,K}"),
        (stages.KΣ, "G^K=-iFA"),
        (stages.I, "I_{\\mathrm{coll}}"),
        (stages.SD, "G^R=D-\\frac{i}{2}A"),
        (stages.R, "\\operatorname{PV}"),
        (stages.N, "F=1+2n"),
        (stages.Q, "q_i\\sim \\pm q_{\\pi(i)}"),
        (stages.C, "C_n(k)="),
    )

    for (value, expected) in displays
        rendered = sprint(show, MIME"text/latex"(), value)
        @test startswith(rendered, "\\[")
        @test occursin(expected, rendered)
        @test endswith(rendered, "\\]")
    end

    @test occursin("Kadanoff-Baym", sprint(show, MIME"text/plain"(), stages.I))
    @test occursin("causal blockers: 0", sprint(show, MIME"text/plain"(), stages.R))
    @test occursin("canonical sectors:", sprint(show, MIME"text/plain"(), stages.C))
end
