using KeldyshContraction, Test
using KeldyshContraction: In, Out, Diagram, Diagrams, Edge
using KeldyshContraction: is_physical, is_conserved, _wick_contraction
using KeldyshContraction: Regularisation.Plus as Plus
using KeldyshContraction: Regularisation.Minus as Minus
import KeldyshContraction as KC
@qfields c::Destroy(Classical) q::Destroy(Quantum)

L_int =
    im * (
        0.5 * c' * q' * (c(Minus) * c(Minus) + q(Minus) * q(Minus)) -
        0.5 * c(Plus) * q(Plus) * (c' * c' + q' * q') +
        c' * q' * (c(Plus) * q(Plus) + c(Minus) * q(Minus))
    )

@testset "vacuum bubble" begin
    @test !iszero(_wick_contraction(L_int; regularise=false))
    @test iszero(_wick_contraction(L_int; regularise=true))
end

@testset "wick contractions first order" begin
    @testset "keldysh Green's function" begin
        using KeldyshContraction: _wick_contraction, regular, In, Out, Diagram
        expr = c(Out()) * c'(In()) * L_int

        @test is_conserved(expr)
        @test is_physical(expr)

        wick_contractions = _wick_contraction(expr.arguments[1].args_nc; regularise=false)
        @test length(wick_contractions) == 4
        regularized_wick = filter(wick_contractions) do cs
            all(regular(c) for c in cs)
        end
        @test length(regularized_wick) == 2
        @test length(unique(Diagram.(regularized_wick))) == 1

        @test isequal(
            Diagram([Edge(c(Out()), q'), Edge(c, c'), Edge(c, c(In())')]),
            set_reg_to_zero(first(keys(_wick_contraction(expr.arguments[1]).diagrams))),
        )

        # ∨ I check these by hand
        simplify = false
        # i (c*c⁻*c⁻*̄c*̄q*̄c)/2
        truth = Diagrams(
            Dict(Diagram([(c(Out()), q'), (c, c'), (c, c'(In()))]) => complex(1.0))
        )
        @test isequal(
            set_reg_to_zero(_wick_contraction(expr.arguments[1]; simplify)), truth
        )

        # i (c*q⁻*q⁻*̄c*̄q*̄c)/2
        truth = Diagrams(
            Dict(Diagram([(c(Out()), q'), (q, c'), (q, c'(In()))]) => complex(1.0))
        )
        @test isequal(
            set_reg_to_zero(_wick_contraction(expr.arguments[2]; simplify)), truth
        )

        # - i( c*c⁺*q⁺*̄c*̄c*̄c) /2
        @test repr(set_reg_to_zero(_wick_contraction(expr.arguments[3]; simplify))) ==
            "-1.0*Gᴷ(x₁,y₁)*Gᴷ(y₁,y₁)*Gᴬ(y₁,x₂)"

        # - i(c*c⁺*q⁺*̄q*̄q*̄ϕ)/2
        @test repr(set_reg_to_zero(_wick_contraction(expr.arguments[4]; simplify))) ==
            "-1.0*Gᴿ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴬ(y₁,x₂)"

        # c*c⁺*q⁺*̄c*̄q*̄c
        truth = Diagrams(
            Dict(
                Diagram([(c(Out()), c'), (c, q'), (q, c'(In()))]) => complex(1.0),
                Diagram([(c(Out()), q'), (c, c'), (q, c'(In()))]) => complex(1.0),
            ),
        )
        @test isequal(
            set_reg_to_zero(_wick_contraction(expr.arguments[5]; simplify)), truth
        )

        # c*c⁻*q⁻*̄c*̄q*̄c
        truth = Diagrams(
            Dict(
                Diagram([(c(Out()), q'), (q, c'), (c, c'(In()))]) => complex(1.0),
                Diagram([(c(Out()), q'), (c, c'), (q, c'(In()))]) => complex(1.0),
            ),
        )
        @test isequal(
            set_reg_to_zero(_wick_contraction(expr.arguments[6]; simplify)), truth
        )

        result = _wick_contraction.(expr.arguments; simplify, _set_reg_to_zero=true)

        diagrams_result = result[1]
        for idx in 2:length(result)
            for (diagram, prefactor) in result[idx]
                push!(diagrams_result, diagram, prefactor)
            end
        end
        diagrams_result

        L = InteractionLagrangian(L_int)
        GF = DressedPropagator(L; simplify=false, _set_reg_to_zero=true)
        @test isequal(GF.keldysh, diagrams_result)
    end

    @testset "quantum-quantum Green's function" begin
        expr = q(Out()) * q'(In()) * L_int

        @test is_conserved(expr)
        @test is_physical(expr)

        # ∨ should this be zero?
        @test !iszero(_wick_contraction(expr; regularise=false))
        @test iszero(_wick_contraction(expr; regularise=true))
    end

    @testset "R/A Green's function first order" begin
        L = InteractionLagrangian(L_int)
        GF = DressedPropagator(L; simplify=false, _set_reg_to_zero=true)

        truth_retarded = Diagrams(
            Dict(
                Diagram([(c(Out()), q'), (q, c'), (c, q'(In()))]) => complex(1.0),
                Diagram([(c(Out()), q'), (c, c'), (c, q'(In()))]) => complex(1.0),
            ),
        )
        truth_advanced = Diagrams(
            Dict(
                Diagram([(q(Out()), c'), (c, q'), (q, c'(In()))]) => complex(1.0),
                Diagram([(q(Out()), c'), (c, c'), (q, c'(In()))]) => complex(-1.0),
            ),
        )
        @test isequal(GF.retarded, truth_retarded)
        @test isequal(GF.advanced, truth_advanced)
    end

    @testset "simplification" begin
        L = InteractionLagrangian(L_int)
        GF_not_simplified = DressedPropagator(L; simplify=false, _set_reg_to_zero=true)
        GF_simplified = DressedPropagator(L; simplify=true, _set_reg_to_zero=true)
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
            GF = DressedPropagator(L; simplify=false, _set_reg_to_zero=true)
            Σ = SelfEnergy(GF)

            kp = Diagram([(c, c')])
            rp = Diagram([(c, q')])
            ap = Diagram([(q, c')])
            advanced_truth = Diagrams(Dict(kp => complex(-1.0), rp => complex(1.0)))
            retarded_truth = Diagrams(Dict(kp => complex(1.0), ap => complex(1.0)))
            keldysh_truth = Diagrams(
                Dict(kp => complex(2.0), rp => complex(-1.0), ap => complex(1.0))
            )

            @test isequal(Σ.advanced, advanced_truth)
            @test isequal(Σ.retarded, retarded_truth)
            @test isequal(Σ.keldysh, keldysh_truth)
            # ^ pretty sure Gerbino et al https://arxiv.org/pdf/2406.20028
            # is wrong and switshes retarded and advanced
            # and we compute the correct with a overall minus sign

            @test isequal(adjoint(Σ.advanced), Σ.retarded)
            @test isequal(adjoint(Σ.keldysh), -1 * Σ.keldysh)

            @testset "simplified" begin
                GF = DressedPropagator(L; simplify=true, _set_reg_to_zero=true)
                Σ = SelfEnergy(GF)
                keldysh_truth = Diagrams(Dict(kp => complex(2.0), rp => complex(-2.0)))
                advanced_truth = Diagrams(Dict(kp => complex(-1.0), rp => complex(1.0)))
                retarded_truth = Diagrams(Dict(kp => complex(1.0), rp => complex(-1.0)))

                @test isequal(Σ.advanced, advanced_truth)
                @test isequal(Σ.retarded, retarded_truth)
                @test isequal(Σ.keldysh, keldysh_truth)
            end
        end
    end

    @testset "Keldysh GF is enough" begin
        using SmallCollections
        using KeldyshContraction: construct_self_energy!, PropagatorType, Diagrams

        L = InteractionLagrangian(L_int)
        GF = DressedPropagator(L; simplify=false, _set_reg_to_zero=true)
        Σ = SelfEnergy(GF)

        expr_K = c(Out()) * c'(In()) * L_int
        G_K1 = _wick_contraction(expr_K; simplify=false, _set_reg_to_zero=true)

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
    GF = DressedPropagator(L, 2, _set_reg_to_zero=true, simplify=true)

    @testset " vacuum" begin
        using KeldyshContraction: filter_nonzero!
        L = InteractionLagrangian(L_int)
        L1 = L(1)
        L2 = L(2)
        vacuum = L1.lagrangian * L2.lagrangian
        expr = _wick_contraction(vacuum; simplify=true, _set_reg_to_zero=true)
        filter_nonzero!(expr)
        @test iszero(expr)
    end

    Σ = SelfEnergy(GF, 2)

    @test_broken isequal(adjoint(Σ.advanced), Σ.retarded)
    @test_broken isequal(adjoint(Σ.keldysh), -1 * Σ.keldysh)

    @test length(topologies(Σ.retarded)[[3]]) ==
        length(topologies(adjoint(Σ.advanced))[[3]])
    @test length(topologies(Σ.retarded)[[3]]) == 4
    @test length(topologies(Σ.retarded)[[2]]) == 6

    @test length(topologies(Σ.keldysh)[[2]]) == 6
    @test length(topologies(Σ.keldysh)[[3]]) == 7
end
