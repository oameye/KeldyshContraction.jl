using KeldyshContraction, Test
import KeldyshContraction as KC

function has_no_zero_loops(component)
    return all(component) do pair
        diagram = first(pair)
        vs = KeldyshContraction.Contraction{Boson}[
            (edge.out, edge.in) for edge in KeldyshContraction.contractions(diagram)
        ]
        return !KeldyshContraction.has_zero_loop(vs)
    end
end

function has_valid_self_energy_diagrams(component)
    return all(component) do pair
        edges = KeldyshContraction.contractions(first(pair))
        return all(KeldyshContraction.is_bulk, edges) &&
               KeldyshContraction.is_irreducible(edges)
    end
end

function irreducible_dressed_component(component)
    out = typeof(component)()
    for (diagram, coefficient) in component
        KC.is_irreducible(KC.contractions(diagram)) || continue
        push!(out, diagram, coefficient)
    end
    return out
end

@testset "number of topologies" begin
    @qfields ϕ::Boson
    c, q = ϕ[Classical], ϕ[Quantum]
    elasctic2boson = -(
        1//2 * (c^2 + q^2) * bar(c) * bar(q) + 1//2 * c * q * (bar(c)^2 + bar(q)^2)
    )
    L_int = InteractionLagrangian(elasctic2boson)

    GF1 = DressedPropagator(L_int, Val(1), Val(3))
    @test length(keys(topologies(GF1.keldysh))) == 1

    GF2 = DressedPropagator(L_int, Val(2), Val(5))
    @test length(keys(topologies(GF2.keldysh))) == 3

    GF3 = DressedPropagator(L_int, Val(3), Val(7))
    @test length(keys(topologies(GF3.keldysh))) == 11
    @test length(unique(sort.(keys(topologies(GF3.keldysh))))) == 8

    irreduciable_topology = []
    for (key, value) in topologies(GF3.keldysh)
        cc = KeldyshContraction.contractions(first(value))
        if KeldyshContraction.is_irreducible(cc)
            push!(irreduciable_topology, key)
        end
    end
    @test length(unique(sort.(irreduciable_topology))) == 5

    @testset "zero-loop filtering" begin
        @test all(has_no_zero_loops, (GF3.keldysh, GF3.retarded, GF3.advanced))
    end

    @testset "third-order static collision path" begin
        Σ = SelfEnergy(GF3)
        @test Σ isa SelfEnergy
        @test all(component -> !isempty(component), (Σ.keldysh, Σ.retarded, Σ.advanced))
        @test all(has_valid_self_energy_diagrams, (Σ.keldysh, Σ.retarded, Σ.advanced))

        # Reuse the already generated order-three object, but Fourier-transform only graphs
        # that can contribute to the self-energy. This avoids a second expensive GF3
        # generation and keeps the regression focused on the downstream static pipeline.
        GF3_irreducible = DressedPropagator(
            irreducible_dressed_component(GF3.keldysh),
            irreducible_dressed_component(GF3.retarded),
            irreducible_dressed_component(GF3.advanced),
            Val(3),
            parameters(GF3),
            target_family(GF3),
        )
        ΣF = SelfEnergy(fourier_transform(GF3_irreducible))
        ΣW = wigner_transform(ΣF; gradient_order=Val(0))
        kinetic = kinetic_expression(ΣW)
        collision = off_shell_collision_expression(kinetic)
        spectral = spectral_dispersive_collision(collision)

        @test KC.order(collision) == 3
        @test KC.order(spectral) == 3
        @test target_family(collision) === ϕ
        @test target_family(spectral) === ϕ

        three_loop_terms = 0
        full_rank_reduced = 0
        partial_rank_reduced = 0
        dependent_three_loop_terms = 0
        dependent_geometries = Dict{Tuple{Int,Int},Int}()
        for expression in
            (collision_offset(spectral), collision_distribution_coefficient(spectral))
            for (term, _) in expression
                loop_frequency_count(term) == 3 || continue
                three_loop_terms += 1

                analysis = @inferred KC.analyze_spectral_dependencies(term, ϕ)
                if KC.has_dependent_shell_support(analysis)
                    dependent_three_loop_terms += 1
                    key = (KC.constraint_rank(analysis), KC.constraint_count(analysis))
                    dependent_geometries[key] = get(dependent_geometries, key, 0) + 1
                end

                rank = spectral_frequency_rank(term)
                reduction = try
                    @inferred KC.general_frequency_reduction(term, ϕ)
                catch error
                    error isa ArgumentError || rethrow()
                    continue
                end
                isempty(reduction) && continue
                if rank == 3
                    full_rank_reduced += 1
                else
                    @test rank < 3
                    partial_rank_reduced += 1
                end
            end
        end

        @test three_loop_terms == 242
        @test full_rank_reduced > 0
        @test partial_rank_reduced > 0
        @test dependent_three_loop_terms == 25
        @test length(dependent_geometries) == 3
        @test get(dependent_geometries, (3, 4), 0) == 7
    end

    GF4 = DressedPropagator(L_int, Val(4), Val(9))
    @test length(keys(topologies(GF4.keldysh))) == 59
    @test length(unique(sort.(keys(topologies(GF4.keldysh))))) == 17

    irreduciable_topology = []
    for (key, value) in topologies(GF4.keldysh)
        cc = KeldyshContraction.contractions(first(value))
        if KeldyshContraction.is_irreducible(cc)
            push!(irreduciable_topology, key)
        end
    end
    @test length(unique(sort.(irreduciable_topology))) == 11

    @testset "zero-loop filtering" begin
        @test all(has_no_zero_loops, (GF4.keldysh, GF4.retarded, GF4.advanced))
    end
end
