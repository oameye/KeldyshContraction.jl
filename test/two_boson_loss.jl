using KeldyshContraction, Test
using KeldyshContraction: In, Out, Diagram, Diagrams, Edge
using KeldyshContraction: set_reg_to_zero
using KeldyshContraction: is_physical, is_conserved, _wick_contraction
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus
import KeldyshContraction as KC
@qfields c::Boson(Classical) q::Boson(Quantum)

L_int =
    im * (
        0.5 * bar(c) * bar(q) * (c(Minus) * c(Minus) + q(Minus) * q(Minus)) -
        0.5 * c(Plus) * q(Plus) * (bar(c) * bar(c) + bar(q) * bar(q)) +
        bar(c) * bar(q) * (c(Plus) * q(Plus) + c(Minus) * q(Minus))
    )

@testset "vacuum bubble" begin
    @test !iszero(_wick_contraction(L_int, Val(2); regularise=false))
    @test iszero(_wick_contraction(L_int, Val(2); regularise=true))
end

@testset "wick contractions first order" begin
    @testset "keldysh Green's function" begin
        using KeldyshContraction: _wick_contraction, regular, In, Out, Diagram
        expr = c(Out()) * bar(c)(In()) * L_int

        @test is_conserved(expr)
        @test is_physical(expr)

        wick_contractions = _wick_contraction(expr.arguments[1].args_nc; regularise=false)
        @test length(wick_contractions) == 4
        regularized_wick = filter(wick_contractions) do cs
            all(regular(cn) for cn in cs)
        end
        @test length(regularized_wick) == 2
        @test length(unique(map(cn -> Diagram(cn, Val(3), Val(0)), regularized_wick))) == 1

        @test isequal(
            Diagram(
                [Edge(c(Out()), bar(q)), Edge(c, bar(c)), Edge(c, bar(c)(In()))],
                Val(3),
                Val(0),
            ),
            set_reg_to_zero(
                first(keys(_wick_contraction(expr.arguments[1], Val(3)).diagrams))
            ),
        )

        simplify = false
        truth = Diagrams(
            Dict(
                Diagram(
                    [(c(Out()), bar(q)), (c, bar(c)), (c, bar(c)(In()))], Val(3), Val(0)
                ) => Complex{Rational{Int}}(1.0),
            ),
        )
        @test isequal(
            set_reg_to_zero(_wick_contraction(expr.arguments[1], Val(3); simplify)), truth
        )

        truth = Diagrams(
            Dict(
                Diagram(
                    [(c(Out()), bar(q)), (q, bar(c)), (q, bar(c)(In()))], Val(3), Val(0)
                ) => Complex{Rational{Int}}(1.0),
            ),
        )
        @test isequal(
            set_reg_to_zero(_wick_contraction(expr.arguments[2], Val(3); simplify)), truth
        )

        @test repr(
            set_reg_to_zero(_wick_contraction(expr.arguments[3], Val(3); simplify))
        ) == "-1//1*Gᴷ(x₁,y₁)*Gᴷ(y₁,y₁)*Gᴬ(y₁,x₂)"

        @test repr(
            set_reg_to_zero(_wick_contraction(expr.arguments[4], Val(3); simplify))
        ) == "-1//1*Gᴿ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴬ(y₁,x₂)"

        truth = Diagrams(
            Dict(
                Diagram(
                    [(c(Out()), bar(c)), (c, bar(q)), (q, bar(c)(In()))], Val(3), Val(0)
                ) => Complex{Rational{Int}}(1.0),
                Diagram(
                    [(c(Out()), bar(q)), (c, bar(c)), (q, bar(c)(In()))], Val(3), Val(0)
                ) => Complex{Rational{Int}}(1.0),
            ),
        )
        @test isequal(
            set_reg_to_zero(_wick_contraction(expr.arguments[5], Val(3); simplify)), truth
        )

        truth = Diagrams(
            Dict(
                Diagram(
                    [(c(Out()), bar(q)), (q, bar(c)), (c, bar(c)(In()))], Val(3), Val(0)
                ) => Complex{Rational{Int}}(1.0),
                Diagram(
                    [(c(Out()), bar(q)), (c, bar(c)), (q, bar(c)(In()))], Val(3), Val(0)
                ) => Complex{Rational{Int}}(1.0),
            ),
        )
        @test isequal(
            set_reg_to_zero(_wick_contraction(expr.arguments[6], Val(3); simplify)), truth
        )

        result = _wick_contraction.(expr.arguments, Val(3); simplify, _set_reg_to_zero=true)

        diagrams_result = result[1]
        for idx in 2:length(result)
            for (diagram, prefactor) in result[idx]
                push!(diagrams_result, diagram, prefactor)
            end
        end

        L = InteractionLagrangian(L_int)
        GF = DressedPropagator(L, Val(1), Val(3); simplify=false, _set_reg_to_zero=true)
        @test isequal(GF.keldysh, diagrams_result)
    end

    @testset "quantum-quantum Green's function" begin
        expr = q(Out()) * bar(q)(In()) * L_int

        @test is_conserved(expr)
        @test is_physical(expr)

        @test !iszero(_wick_contraction(expr, Val(3); regularise=false))
        @test iszero(_wick_contraction(expr, Val(3); regularise=true))
    end

    @testset "R/A Green's function first order" begin
        L = InteractionLagrangian(L_int)
        GF = DressedPropagator(L, Val(1), Val(3); simplify=false, _set_reg_to_zero=true)

        truth_retarded = Diagrams(
            Dict(
                Diagram(
                    [(c(Out()), bar(q)), (q, bar(c)), (c, bar(q)(In()))], Val(3), Val(0)
                ) => Complex{Rational{Int}}(1.0),
                Diagram(
                    [(c(Out()), bar(q)), (c, bar(c)), (c, bar(q)(In()))], Val(3), Val(0)
                ) => Complex{Rational{Int}}(1.0),
            ),
        )
        truth_advanced = Diagrams(
            Dict(
                Diagram(
                    [(q(Out()), bar(c)), (c, bar(q)), (q, bar(c)(In()))], Val(3), Val(0)
                ) => Complex{Rational{Int}}(1.0),
                Diagram(
                    [(q(Out()), bar(c)), (c, bar(c)), (q, bar(c)(In()))], Val(3), Val(0)
                ) => Complex{Rational{Int}}(-1.0),
            ),
        )
        @test isequal(set_reg_to_zero(GF.retarded), truth_retarded)
        @test isequal(set_reg_to_zero(GF.advanced), truth_advanced)
    end

    @testset "simplification" begin
        L = InteractionLagrangian(L_int)
        GF_not_simplified = DressedPropagator(
            L, Val(1), Val(3); simplify=false, _set_reg_to_zero=true
        )
        GF_simplified = DressedPropagator(
            L, Val(1), Val(3); simplify=true, _set_reg_to_zero=true
        )
        collect(keys(GF_simplified.keldysh.diagrams))
        collect(values(GF_simplified.keldysh.diagrams))
        collect(keys(GF_not_simplified.keldysh.diagrams))
        collect(values(GF_not_simplified.keldysh.diagrams))
    end
end

@testset "self-energy first order" begin
    using KeldyshContraction: Edge, matrix, Diagrams, Diagram
    L = InteractionLagrangian(L_int)

    @testset "correctness check" begin
        @testset "first order" begin
            GF = DressedPropagator(L, Val(1), Val(3); simplify=false, _set_reg_to_zero=true)
            Σ = SelfEnergy(GF, Val(1))

            kp = Diagram([(c, bar(c))], Val(1), Val(0))
            rp = Diagram([(c, bar(q))], Val(1), Val(0))
            ap = Diagram([(q, bar(c))], Val(1), Val(0))
            advanced_truth = Diagrams(
                Dict(kp => Complex{Rational{Int}}(-1.0), rp => Complex{Rational{Int}}(1.0))
            )
            retarded_truth = Diagrams(
                Dict(kp => Complex{Rational{Int}}(1.0), ap => Complex{Rational{Int}}(1.0))
            )
            keldysh_truth = Diagrams(
                Dict(
                    kp => Complex{Rational{Int}}(2.0),
                    rp => Complex{Rational{Int}}(-1.0),
                    ap => Complex{Rational{Int}}(1.0),
                ),
            )

            @test isequal(set_reg_to_zero(Σ.advanced), advanced_truth)
            @test isequal(set_reg_to_zero(Σ.retarded), retarded_truth)
            @test isequal(set_reg_to_zero(Σ.keldysh), keldysh_truth)

            @test isequal(adjoint(set_reg_to_zero(Σ.advanced)), set_reg_to_zero(Σ.retarded))
            @test isequal(
                adjoint(set_reg_to_zero(Σ.keldysh)), -1 * set_reg_to_zero(Σ.keldysh)
            )

            @testset "simplified" begin
                GF = DressedPropagator(
                    L, Val(1), Val(3); simplify=true, _set_reg_to_zero=true
                )
                Σ = SelfEnergy(GF, Val(1))
                keldysh_truth = Diagrams(
                    Dict(
                        kp => Complex{Rational{Int}}(2.0),
                        rp => Complex{Rational{Int}}(-2.0),
                    ),
                )
                advanced_truth = Diagrams(
                    Dict(
                        kp => Complex{Rational{Int}}(-1.0),
                        rp => Complex{Rational{Int}}(1.0),
                    ),
                )
                retarded_truth = Diagrams(
                    Dict(
                        kp => Complex{Rational{Int}}(1.0),
                        rp => Complex{Rational{Int}}(-1.0),
                    ),
                )

                @test isequal(set_reg_to_zero(Σ.advanced), advanced_truth)
                @test isequal(set_reg_to_zero(Σ.retarded), retarded_truth)
                @test isequal(set_reg_to_zero(Σ.keldysh), keldysh_truth)
            end
        end
    end

    @testset "Keldysh GF is enough" begin
        using SmallCollections
        using KeldyshContraction: construct_self_energy!, PropagatorType, Diagrams

        L = InteractionLagrangian(L_int)
        GF = DressedPropagator(L, Val(1), Val(3); simplify=false, _set_reg_to_zero=true)
        Σ = SelfEnergy(GF, Val(1))

        expr_K = c(Out()) * bar(c)(In()) * L_int
        G_K1 = _wick_contraction(expr_K, Val(3); simplify=false, _set_reg_to_zero=true)

        self_energy = SmallCollections.SmallDict{3,PropagatorType.T,Diagrams}((
            PropagatorType.Advanced => Diagrams{1,0}(),
            PropagatorType.Retarded => Diagrams{1,0}(),
            PropagatorType.Keldysh => Diagrams{1,0}(),
        ))
        construct_self_energy!(self_energy, G_K1)
        @test isequal(self_energy[PropagatorType.Advanced], Σ.advanced)
        @test isequal(self_energy[PropagatorType.Retarded], Σ.retarded)
    end
end

@testset "second order" begin
    L = InteractionLagrangian(L_int)
    GF = DressedPropagator(L, Val(2), Val(5), _set_reg_to_zero=true, simplify=true)

    @testset "vacuum" begin
        using KeldyshContraction: filter_nonzero!
        L = InteractionLagrangian(L_int)
        L1 = L(1)
        L2 = L(2)
        vacuum = L1.lagrangian * L2.lagrangian
        expr = _wick_contraction(vacuum, Val(4); simplify=true, _set_reg_to_zero=true)
        filter_nonzero!(expr)
        @test iszero(expr)
    end

    Σ = SelfEnergy(GF, Val(2))

    # 9 of the 11 diagrams line up. The two that do not differ only in which leg carries
    # the regularisation of an equal-time tadpole: Gᴿ(y⁺,y) against Gᴿ(y,y⁻), which have
    # the same `subtraction` and so are the same propagator written two ways. Wick
    # contraction emits both spellings, so closing this needs a canonical form for the
    # regularisation of an equal-position edge, not a change to the adjoint.
    @test_broken isequal(adjoint(Σ.advanced), Σ.retarded)
    @test_broken isequal(adjoint(Σ.keldysh), -1 * Σ.keldysh)

    @test length(topologies(Σ.retarded)[[3]]) ==
        length(topologies(adjoint(Σ.advanced))[[3]])
    @test length(topologies(Σ.retarded)[[3]]) == 4
    @test length(topologies(Σ.retarded)[[2]]) == 7

    @test length(topologies(Σ.keldysh)[[2]]) == 7
    @test length(topologies(Σ.keldysh)[[3]]) == 7
end
