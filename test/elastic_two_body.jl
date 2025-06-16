using KeldyshContraction, Test
using KeldyshContraction: In, Out, _wick_contraction
using OrderedCollections
using Random
Random.seed!(1234) # for reproducibility

@qfields c::Destroy(Classical) q::Destroy(Quantum)
elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
L_int = InteractionLagrangian(elasctic2boson)

@testset "first order" begin
    @testset "Bubble diagrams" begin
        using KeldyshContraction: filter_nonzero!
        expr = _wick_contraction(elasctic2boson; simplify=true)
        filter_nonzero!(expr)
        @test iszero(expr)
    end

    @testset "wick contraction" begin
        using KeldyshContraction: _wick_contraction, regular, In, Out, Diagram, Diagrams

        expr = c(Out()) * q'(In()) * elasctic2boson

        @test KeldyshContraction.is_conserved(expr)
        @test KeldyshContraction.is_physical(expr)

        # ∨ I check these by hand
        # 0.5*(c*c*c*̄q*̄c*̄q
        truth = Diagrams(
            Dict(
                Diagram([(c(Out()), c'), (c, q'), (c, q(In())')]) => -0.0 + 1.0 * im,
                Diagram([(c(Out()), q'), (c, c'), (c, q(In())')]) => -0.0 + 1.0 * im,
            ),
        )
        result = _wick_contraction(expr.arguments[1])
        KeldyshContraction._simplify_prefactors!(result)
        @test isequal(result, truth)
        # The keldysh in and keldysh out will disappear later
    end

    @testset "green's function" begin
        # using KeldyshContraction: _conj
        L = InteractionLagrangian(elasctic2boson)
        GF = DressedPropagator(L)

        @test_broken iszero(_conj(GF.advanced) - GF.retarded)
        @test_broken isequal(KeldyshContraction._conj(GF.keldysh), -1 * GF.keldysh)
        # TODO this is not equal to do G^R and G^A being at different coordinates
        # switch in and out coordinate when adjoint(::DressedPropagator)
    end

    @testset "self-energy" begin
        # using KeldyshContraction: _conj
        L = InteractionLagrangian(elasctic2boson)
        GF = DressedPropagator(L)
        Σ = SelfEnergy(GF)

        @test iszero(Σ.keldysh)
        @test isequal(adjoint(Σ.advanced), Σ.retarded)
        @test isequal(adjoint(Σ.keldysh), -1 * Σ.keldysh)
    end
end

@testset "second order" begin
    @testset "vacuum" begin
        L1 = L_int(1)
        L2 = L_int(2)

        using KeldyshContraction: filter_nonzero!

        vacuum = L1.lagrangian * L2.lagrangian
        map(vacuum.arguments) do arg
            _wick_contraction(arg; simplify=false)
        end
        _wick_contraction(vacuum.arguments[1]; simplify=false)

        expr = _wick_contraction(vacuum; simplify=true)
        filter_nonzero!(expr)
        @test iszero(expr)
    end

    @testset "individual terms" begin
        using KeldyshContraction: Diagram, Diagrams, Bulk, Edge
        import KeldyshContraction as KC
        using Combinatorics
        order = 2
        in_out = c(Out()) * c'(In())
        l = length(L_int.lagrangian)

        E = KC.number_of_propagators(L_int) * order + 1 # +1 for in_out
        dict = OrderedDict()
        keys = []
        prefactor = -1 * im * im^order / factorial(order)

        for coefficients in Combinatorics.multiexponents(l, order)
            idxs = KC.indices_from_counts(coefficients) # will be of length order
            mult = Combinatorics.multinomial(coefficients...)
            qmul = mult * prod(L_int(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))

            term = prefactor * in_out * qmul
            diagrams = _wick_contraction(term; simplify=true)
            KC._simplify_prefactors!(diagrams)
            dict[term] = diagrams
            push!(keys, term)
        end
        dict

        @test length(dict) == 10
        for idx in [5, 7, 10]
            @test isempty(dict[keys[idx]].diagrams)
        end

        truth1 = Diagrams(
            Dict(
                # -1.0*Gᴿ(x₁,y₁)*Gᴷ(y₁,y₁)*Gᴿ(y₁,y₂)*Gᴷ(y₂,y₂)*Gᴷ(y₂,x₂)
                Diagram([
                    (c(Out()), q'(Bulk(1))),     # Gᴿ(x₁,y₁)
                    (c(Bulk(1)), c'(Bulk(1))),   # Gᴷ(y₁,y₁)
                    (c(Bulk(1)), q'(Bulk(2))),   # Gᴿ(y₁,y₂)
                    (c(Bulk(2)), c'(Bulk(2))),   # Gᴷ(y₂,y₂)
                    (c(Bulk(2)), c'(In())),      # Gᴷ(y₂,x₂)
                ]) => -1.0 + 0im,
                # -1.0*Gᴷ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴿ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴷ(y₁,x₂)
                Diagram([
                    (c(Out()), c'(Bulk(1))),     # Gᴷ(x₁,y₁)
                    (c(Bulk(1)), c'(Bulk(2))),   # Gᴷ(y₁,y₂)
                    (c(Bulk(2)), q'(Bulk(1))),   # Gᴿ(y₂,y₁)
                    (c(Bulk(2)), q'(Bulk(2))),   # Gᴿ(y₂,y₂)
                    (c(Bulk(1)), c'(In())),      # Gᴷ(y₁,x₂)
                ]) => -1.0 + 0im,
                # -1.0*Gᴿ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴿ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴷ(y₂,x₂)
                Diagram([
                    (c(Out()), q'(Bulk(1))),     # Gᴿ(x₁,y₁)
                    (c(Bulk(1)), c'(Bulk(2))),   # Gᴷ(y₁,y₂)
                    (c(Bulk(1)), q'(Bulk(2))),   # Gᴿ(y₁,y₂)
                    (c(Bulk(2)), c'(Bulk(1))),   # Gᴷ(y₂,y₁)
                    (c(Bulk(2)), c'(In())),      # Gᴷ(y₂,x₂)
                ]) => -1.0 + 0im,
                # -1.0*Gᴿ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴷ(y₁,x₂)
                Diagram([
                    (c(Out()), q'(Bulk(1))),     # Gᴿ(x₁,y₁)
                    (c(Bulk(1)), c'(Bulk(2))),   # Gᴷ(y₁,y₂)
                    (c(Bulk(2)), c'(Bulk(1))),   # Gᴷ(y₂,y₁)
                    (c(Bulk(2)), q'(Bulk(2))),   # Gᴿ(y₂,y₂)
                    (c(Bulk(1)), c'(In())),      # Gᴷ(y₁,x₂)
                ]) => -1.0 + 0im,
                # -1.0*Gᴿ(x₁,y₁)*Gᴷ(y₁,y₁)*Gᴷ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴷ(y₂,x₂)
                Diagram([
                    (c(Out()), q'(Bulk(1))),     # Gᴿ(x₁,y₁)
                    (c(Bulk(1)), c'(Bulk(1))),   # Gᴷ(y₁,y₁)
                    (c(Bulk(1)), c'(Bulk(2))),   # Gᴷ(y₁,y₂)
                    (c(Bulk(2)), q'(Bulk(2))),   # Gᴿ(y₂,y₂)
                    (c(Bulk(2)), c'(In())),      # Gᴷ(y₂,x₂)
                ]) => -1.0 + 0im,
                # -1.0*Gᴿ(x₁,y₁)*Gᴿ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴷ(y₂,y₂)*Gᴷ(y₁,x₂)
                Diagram([
                    (c(Out()), q'(Bulk(1))),     # Gᴿ(x₁,y₁)
                    (c(Bulk(1)), q'(Bulk(2))),   # Gᴿ(y₁,y₂)
                    (c(Bulk(2)), c'(Bulk(1))),   # Gᴷ(y₂,y₁)
                    (c(Bulk(2)), c'(Bulk(2))),   # Gᴷ(y₂,y₂)
                    (c(Bulk(1)), c'(In())),      # Gᴷ(y₁,x₂)
                ]) => -1.0 + 0im,
                # -1.0*Gᴷ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴷ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴷ(y₂,x₂)
                Diagram([
                    (c(Out()), c'(Bulk(1))),     # Gᴷ(x₁,y₁)
                    (c(Bulk(1)), q'(Bulk(1))),   # Gᴿ(y₁,y₁)
                    (c(Bulk(1)), c'(Bulk(2))),   # Gᴷ(y₁,y₂)
                    (c(Bulk(2)), q'(Bulk(2))),   # Gᴿ(y₂,y₂)
                    (c(Bulk(2)), c'(In())),      # Gᴷ(y₂,x₂)
                ]) => -1.0 + 0im,
                # -1.0*Gᴷ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴿ(y₁,y₂)*Gᴷ(y₂,y₂)*Gᴷ(y₂,x₂)
                Diagram([
                    (c(Out()), c'(Bulk(1))),     # Gᴷ(x₁,y₁)
                    (c(Bulk(1)), q'(Bulk(1))),   # Gᴿ(y₁,y₁)
                    (c(Bulk(1)), q'(Bulk(2))),   # Gᴿ(y₁,y₂)
                    (c(Bulk(2)), c'(Bulk(2))),   # Gᴷ(y₂,y₂)
                    (c(Bulk(2)), c'(In())),      # Gᴷ(y₂,x₂)
                ]) => -1.0 + 0im,
            ),
        )
        @test isequal(dict[keys[1]], truth1)

        truth2 = Diagrams(
            Dict(
                # Gᴿ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴬ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴷ(y₂,x₂)
                Diagram([
                    (c(Out()), q'(Bulk(1))),
                    (c(Bulk(1)), q'(Bulk(1))),
                    (q(Bulk(1)), c'(Bulk(2))),
                    (c(Bulk(2)), q'(Bulk(2))),
                    (c(Bulk(2)), c'(In())),
                ]) => 1.0 + 0im,
                # -1.0*Gᴿ(x₁,y₁)*Gᴬ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴬ(y₁,x₂)
                Diagram([
                    (c(Out()), q'(Bulk(1))),
                    (q(Bulk(1)), c'(Bulk(2))),
                    (c(Bulk(2)), c'(Bulk(1))),
                    (c(Bulk(2)), q'(Bulk(2))),
                    (q(Bulk(1)), c'(In())),
                ]) => -1.0 + 0im,
                # -1.0*Gᴿ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴿ(y₁,y₂)*Gᴬ(y₂,y₁)*Gᴬ(y₂,x₂)
                Diagram([
                    (c(Out()), q'(Bulk(1))),
                    (c(Bulk(1)), c'(Bulk(2))),
                    (c(Bulk(1)), q'(Bulk(2))),
                    (q(Bulk(2)), c'(Bulk(1))),
                    (q(Bulk(2)), c'(In())),
                ]) => -1.0 + 0im,
                # Gᴿ(x₁,y₁)*Gᴷ(y₁,y₁)*Gᴿ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴬ(y₂,x₂)
                Diagram([
                    (c(Out()), q'(Bulk(1))),
                    (c(Bulk(1)), c'(Bulk(1))),
                    (c(Bulk(1)), q'(Bulk(2))),
                    (c(Bulk(2)), q'(Bulk(2))),
                    (q(Bulk(2)), c'(In())),
                ]) => 1.0 + 0im,
                # -1.0*Gᴷ(x₁,y₁)*Gᴬ(y₁,y₂)*Gᴿ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴬ(y₁,x₂)
                Diagram([
                    (c(Out()), c'(Bulk(1))),
                    (q(Bulk(1)), c'(Bulk(2))),
                    (c(Bulk(2)), q'(Bulk(1))),
                    (c(Bulk(2)), q'(Bulk(2))),
                    (q(Bulk(1)), c'(In())),
                ]) => -1.0 + 0im,
                # Gᴿ(x₁,y₁)*Gᴿ(y₁,y₂)*Gᴬ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴷ(y₁,x₂)
                Diagram([
                    (c(Out()), q'(Bulk(1))),
                    (c(Bulk(1)), q'(Bulk(2))),
                    (q(Bulk(2)), c'(Bulk(1))),
                    (c(Bulk(2)), q'(Bulk(2))),
                    (c(Bulk(1)), c'(In())),
                ]) => 1.0 + 0im,
                # Gᴷ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴿ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴬ(y₂,x₂)
                Diagram([
                    (c(Out()), c'(Bulk(1))),
                    (c(Bulk(1)), q'(Bulk(1))),
                    (c(Bulk(1)), q'(Bulk(2))),
                    (c(Bulk(2)), q'(Bulk(2))),
                    (q(Bulk(2)), c'(In())),
                ]) => 1.0 + 0im,
            ),
        )
        @test isequal(dict[keys[2]], truth2)

        # @test repr(dict[keys[3]]) ==
        #     "-1.0*Gᴷ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴷ(y₁,y₂)*Gᴷ(y₂,y₂)*Gᴬ(y₂,x₂) + -0.5*Gᴷ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴷ(y₁,y₂)*Gᴿ(y₂,y₁)*Gᴬ(y₂,x₂) + 2.0*Gᴷ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴷ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴷ(y₂,x₂) + Gᴷ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴿ(y₁,y₂)*Gᴷ(y₂,y₂)*Gᴷ(y₂,x₂) + -1.0*Gᴷ(x₁,y₁)*Gᴬ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴷ(y₁,x₂) + Gᴿ(x₁,y₁)*Gᴷ(y₁,y₁)*Gᴷ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴷ(y₂,x₂) + -1.0*Gᴷ(x₁,y₁)*Gᴿ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴷ(y₂,y₂)*Gᴬ(y₁,x₂) + -0.5*Gᴿ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴷ(y₁,y₂)*Gᴬ(y₂,y₁)*Gᴷ(y₂,x₂) + -1.0*Gᴷ(x₁,y₁)*Gᴷ(y₁,y₁)*Gᴬ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴷ(y₂,x₂) + -1.0*Gᴷ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴬ(y₁,x₂) + -1.0*Gᴿ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴬ(y₂,y₁)*Gᴷ(y₂,y₂)*Gᴷ(y₁,x₂) + Gᴷ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴿ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴷ(y₁,x₂) + Gᴿ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴷ(y₁,x₂) + -1.0*Gᴿ(x₁,y₁)*Gᴷ(y₁,y₁)*Gᴷ(y₁,y₂)*Gᴷ(y₂,y₂)*Gᴬ(y₂,x₂) + -0.5*Gᴿ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴷ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴬ(y₂,x₂)"

        # @test repr(dict[keys[4]]) ==
        #     "-1.0*Gᴿ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴬ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴷ(y₂,x₂) + -1.0*Gᴿ(x₁,y₁)*Gᴷ(y₁,y₁)*Gᴿ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴬ(y₂,x₂) + -1.0*Gᴿ(x₁,y₁)*Gᴿ(y₁,y₂)*Gᴬ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴷ(y₁,x₂) + -0.5*Gᴿ(x₁,y₁)*Gᴿ(y₁,y₂)*Gᴿ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴬ(y₂,x₂) + -1.0*Gᴿ(x₁,y₁)*Gᴬ(y₁,y₂)*Gᴿ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴷ(y₁,x₂) + -1.0*Gᴷ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴿ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴬ(y₂,x₂) + -1.0*Gᴿ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴿ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴬ(y₁,x₂) + -0.5*Gᴿ(x₁,y₁)*Gᴿ(y₁,y₂)*Gᴿ(y₁,y₂)*Gᴬ(y₂,y₁)*Gᴷ(y₂,x₂)"

        # @test repr(dict[keys[6]]) ==
        #     "-0.5*Gᴷ(x₁,y₁)*Gᴬ(y₁,y₂)*Gᴬ(y₁,y₂)*Gᴿ(y₂,y₁)*Gᴬ(y₂,x₂) + Gᴿ(x₁,y₁)*Gᴬ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴬ(y₁,x₂) + -1.0*Gᴿ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴬ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴷ(y₂,x₂) + Gᴿ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴬ(y₁,y₂)*Gᴷ(y₂,y₂)*Gᴬ(y₂,x₂) + -0.5*Gᴿ(x₁,y₁)*Gᴬ(y₁,y₂)*Gᴬ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴬ(y₂,x₂) + Gᴷ(x₁,y₁)*Gᴬ(y₁,y₂)*Gᴿ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴬ(y₁,x₂) + Gᴷ(x₁,y₁)*Gᴿ(y₁,y₂)*Gᴬ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴬ(y₁,x₂) + -1.0*Gᴷ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴿ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴬ(y₂,x₂)"

        # @test repr(dict[keys[8]]) ==
        #     "-1.0*Gᴷ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴬ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴬ(y₂,x₂) + Gᴷ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴷ(y₁,y₂)*Gᴷ(y₂,y₂)*Gᴬ(y₂,x₂) + -1.0*Gᴷ(x₁,y₁)*Gᴷ(y₁,y₁)*Gᴬ(y₁,y₂)*Gᴷ(y₂,y₂)*Gᴬ(y₂,x₂) + Gᴷ(x₁,y₁)*Gᴷ(y₁,y₁)*Gᴬ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴷ(y₂,x₂) + Gᴷ(x₁,y₁)*Gᴬ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴷ(y₁,x₂) + -1.0*Gᴷ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴬ(y₂,y₁)*Gᴷ(y₂,y₂)*Gᴬ(y₁,x₂) + -1.0*Gᴷ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴷ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴷ(y₂,x₂) + Gᴷ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴷ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴬ(y₁,x₂)"

        # @test repr(dict[keys[9]]) ==
        #     "Gᴿ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴬ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴷ(y₂,x₂) + -1.0*Gᴿ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴬ(y₁,y₂)*Gᴿ(y₂,y₁)*Gᴬ(y₂,x₂) + -1.0*Gᴿ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴬ(y₁,y₂)*Gᴷ(y₂,y₂)*Gᴬ(y₂,x₂) + -1.0*Gᴷ(x₁,y₁)*Gᴿ(y₁,y₂)*Gᴬ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴬ(y₁,x₂) + Gᴿ(x₁,y₁)*Gᴬ(y₁,y₂)*Gᴿ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴷ(y₁,x₂) + Gᴷ(x₁,y₁)*Gᴿ(y₁,y₁)*Gᴿ(y₁,y₂)*Gᴿ(y₂,y₂)*Gᴬ(y₂,x₂) + Gᴿ(x₁,y₁)*Gᴷ(y₁,y₂)*Gᴿ(y₂,y₁)*Gᴿ(y₂,y₂)*Gᴬ(y₁,x₂)"
    end

    GF = DressedPropagator(L_int; order=2)

    @testset "collective terms" begin
        GF = DressedPropagator(L_int; order=2)
        GF.retarded.diagrams

        terms_k = collect(keys(GF.keldysh.diagrams))
        @test length(terms_k) == 18

        terms_r = collect(GF.retarded.diagrams)
        @test length(terms_r) == 6

        terms_a = collect(GF.advanced.diagrams)
        @test length(terms_a) == 6
    end

    Σ = SelfEnergy(GF; order=2)

    @test_broken isequal(adjoint(Σ.advanced), Σ.retarded) # up to some swap
    @test isequal(adjoint(Σ.keldysh), -1 * Σ.keldysh)
end
