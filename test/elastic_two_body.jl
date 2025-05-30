using KeldyshContraction, Test
using KeldyshContraction: In, Out

@qfields c::Destroy(Classical) q::Destroy(Quantum)
elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
L_int = InteractionLagrangian(elasctic2boson)

# using KeldyshContraction: Minus, Plus
# elasctic2boson_reguralize = 0.5 * (c(Minus)^2 + q(Minus)^2) * c' * q' + 0.5 * c(Plus) * q(Plus) * ((c')^2 + (q')^2)
# L_int = InteractionLagrangian(elasctic2boson_reguralize)

@testset "first order" begin
    @testset "Bubble diagrams" begin
        using KeldyshContraction: filter_nonzero!
        expr = wick_contraction(elasctic2boson; simplify=true)
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
        result = wick_contraction(expr.arguments[1])
        KeldyshContraction._simplify!(result)
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
    L1 = L_int(1)
    L2 = L_int(2)

    @testset "vacuum" begin
        using KeldyshContraction: filter_nonzero!

        vacuum = L1.lagrangian * L2.lagrangian
        map(vacuum.arguments) do arg
            wick_contraction(arg; simplify=false)
        end
        wick_contraction(vacuum.arguments[1]; simplify=false)

        expr = wick_contraction(vacuum; simplify=true)
        filter_nonzero!(expr)
        @test iszero(expr)
    end

    expr = c(Out()) * c'(In()) * L1.lagrangian * L2.lagrangian

    import KeldyshContraction as KC
    using Combinatorics
    order = 2
    in_out = c(Out()) * c'(In())
    l = length(L_int.lagrangian)

    E = KC.number_of_propagators(L_int)*order + 1 # +1 for in_out
    dict = Dict()
    prefactor = -1 * im * im^order / factorial(order)

    for coefficients in Combinatorics.multiexponents(l, order)
        idxs = KC.expand_coefficients(coefficients) # will be of length order
        mult = Combinatorics.multinomial(coefficients...)
        qmul = mult * prod(L_int(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))

        term = prefactor * in_out * qmul
        diagrams = wick_contraction(term; simplify=true)
        dict[term] = diagrams
    end
    dict
    i=5
    @show collect(dict)[i][1]
    collect(dict)[i][2].diagrams

    for coefficients in Combinatorics.multiexponents(l, order)
        mult = Combinatorics.multinomial(coefficients...)
        @show (coefficients, mult)
    end
    # ∨ I check these by hand
    # 0.25*(c*c*c*c*c*̄q*̄c*̄q*̄c*̄c)
    # truth = Diagrams(
    #     Dict(
    #         Diagram([(c(Out()), c'), (c, q'), (c, q'), (c, q'), (c, q(In())')]) => 0.0 - 2.0 * im,
    #         Diagram([(c(Out()), q'), (c, c'), (c, q(In())')]) => 0.0 - 2.0 * im,
    #     ),
    # )
    # wick_contraction(expr.arguments[1]).diagrams
    # ^ TODO we should write this test
    # The keldysh in and keldysh out will disappear later

    GF = DressedPropagator(L_int; order=2)
    GF.retarded.diagrams

    terms_k = collect(keys(GF.keldysh.diagrams))
    @test length(terms_k) == 18

    terms_r = collect(GF.retarded.diagrams)
    @test length(terms_r) == 6

    terms_a = collect(GF.advanced.diagrams)
    @test length(terms_a) == 6

    terms_a[2]
    KeldyshContraction.adjoint_diagram(terms_r[4])

    Σ = SelfEnergy(GF; order=2)

    @test_broken isequal(adjoint(Σ.advanced), Σ.retarded) # up to some swap
    @test isequal(adjoint(Σ.keldysh), -1 * Σ.keldysh)
end
