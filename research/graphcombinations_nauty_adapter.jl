using BenchmarkTools
using KeldyshContraction
using Test

import GraphCombinations as GC
import KeldyshContraction as KC

const GC_SHA = "ec5f84d4b1b1f21b2dacdb00ce4ed9af6bd68fa8"
const Out = KC.Out
const In = KC.In
const Bulk = KC.Bulk
const OccupationPolynomial = KC.OccupationPolynomial
const OccupationAtom = KC.OccupationAtom

function gc_result(edges::Vector{Tuple{Int,Int}}, colors::Vector{Int})
    graph = GC.DirectedGCGraph(
        Pair{Int,Int}[source => target for (source, target) in edges], length(colors)
    )
    return GC.canonicalize_directed(graph, colors)
end

function gc_position_edges(vs, graph_positions)
    edges = Tuple{Int,Int}[]
    sizehint!(edges, length(vs))
    for item in vs
        out, in = KC.positions(item)
        push!(
            edges,
            (
                KC.position_vertex(graph_positions, out),
                KC.position_vertex(graph_positions, in),
            ),
        )
    end
    return edges
end

function gc_topology_result(vs, graph_positions=KC.canonicalization_positions(vs))
    edges = unique(gc_position_edges(vs, graph_positions))
    colors = KC.position_labels(graph_positions)
    return gc_result(edges, colors)
end

function gc_physical_result(vs, graph_positions=KC.canonicalization_positions(vs))
    direct_edges = gc_position_edges(vs, graph_positions)
    simple = length(unique(direct_edges)) == length(direct_edges)
    if simple && KC.uniform_coloring(vs)
        return gc_result(direct_edges, KC.position_labels(graph_positions))
    end

    physical_colors = KC.propagator_colors(vs)
    npositions = length(graph_positions)
    labels = Vector{Int}(undef, npositions + length(vs))
    copyto!(labels, 1, KC.position_labels(graph_positions), 1, npositions)
    edges = Tuple{Int,Int}[]
    sizehint!(edges, 2 * length(vs))
    for (i, item) in enumerate(vs)
        color_index = searchsortedfirst(physical_colors, KC.propagator_color(item))
        edge_vertex = npositions + i
        labels[edge_vertex] = 3 + color_index
        out, in = KC.positions(item)
        push!(edges, (KC.position_vertex(graph_positions, out), edge_vertex))
        push!(edges, (edge_vertex, KC.position_vertex(graph_positions, in)))
    end
    return gc_result(edges, labels)
end

function gc_canonicalize(vs::Vector{T}) where {T}
    isempty(vs) && return copy(vs)
    graph_positions = KC.canonicalization_positions(vs)
    result = gc_physical_result(vs, graph_positions)
    old_to_canonical = GC.vertex_mapping(GC.canonical_relabeling(result))
    canonical_to_old = sortperm(old_to_canonical)
    mapping = KC.make_permutation_dict(canonical_to_old, graph_positions, vs)
    return T[KC.relabel_bulk_positions(item, mapping) for item in vs]
end

function gc_legacy_topology(vs, ::Val{E2}) where {E2}
    isempty(vs) && return KC.bulk_multiplicity(Tuple{Int8,Int8}[], Val(E2))

    position_pairs = Tuple{Int8,Int8}[KC.integer_positions(item) for item in vs]
    flattened = collect(Iterators.flatten(position_pairs))
    max_label = length(unique(flattened))
    has_out = typemin(Int8) in flattened

    edges = Tuple{Int,Int}[]
    for pair in position_pairs
        vertices = if typemin(Int8) in pair
            (1, Int(last(pair)) + Int(has_out))
        elseif typemax(Int8) in pair
            (Int(first(pair)) + Int(has_out), max_label)
        else
            (Int(first(pair)) + Int(has_out), Int(last(pair)) + Int(has_out))
        end
        vertices in edges || push!(edges, vertices)
    end

    result = gc_result(edges, ones(Int, max_label))
    old_to_canonical = GC.vertex_mapping(GC.canonical_relabeling(result))
    canonical_to_old = sortperm(old_to_canonical)
    mapping = KC.legacy_topology_permutation_dict(canonical_to_old, max_label, has_out)
    topology_edges = Tuple{Int8,Int8}[
        KC.integer_positions(KC.relabel_bulk_positions(item, mapping)) for item in vs
    ]
    return KC.bulk_multiplicity(topology_edges, Val(E2))
end

function gc_loop_graph_data(
    sector::KC.ReducedCollisionSector{S}, monomial::KC.OccupationMonomial{S}
) where {S<:KC.Statistics}
    basis = KC.momentum_basis(sector)
    external = KC.external_wigner_momentum(sector)
    external_index = KC._external_basis_index(basis, external)
    loop_indices = KC._loop_basis_indices(basis, external)
    nloops = length(loop_indices)

    builder = KC._LoopCanonicalGraphBuilder()
    root = KC._add_loop_vertex!(builder, KC._loop_graph_color(1))
    pair_vertices = Vector{Int}(undef, nloops)
    positive_vertices = Vector{Int}(undef, nloops)
    negative_vertices = Vector{Int}(undef, nloops)
    loop_incidences = [Tuple{Int,KC.MomentumCoefficient}[] for _ in 1:nloops]

    for slot in 1:nloops
        pair = KC._add_loop_vertex!(builder, KC._loop_graph_color(2))
        positive = KC._add_loop_vertex!(builder, KC._loop_graph_color(3))
        negative = KC._add_loop_vertex!(builder, KC._loop_graph_color(3))
        pair_vertices[slot] = pair
        positive_vertices[slot] = positive
        negative_vertices[slot] = negative
        KC._add_loop_edge!(builder, root, pair)
        KC._add_loop_edge!(builder, pair, positive)
        KC._add_loop_edge!(builder, pair, negative)
    end

    KC._add_loop_semantics!(
        builder,
        root,
        sector,
        monomial,
        external_index,
        loop_indices,
        positive_vertices,
        negative_vertices,
        loop_incidences,
    )

    color_classes = sort!(unique(copy(builder.colors)))
    labels = Int[searchsortedfirst(color_classes, color) for color in builder.colors]
    result = gc_result(builder.edges, labels)
    return result, pair_vertices, loop_incidences, basis, external, nloops
end

function gc_canonical_loop_transform(
    sector::KC.ReducedCollisionSector{S}, monomial::KC.OccupationMonomial{S}
) where {S<:KC.Statistics}
    result, pair_vertices, loop_incidences, basis, external, nloops = gc_loop_graph_data(
        sector, monomial
    )

    rank = GC.vertex_mapping(GC.canonical_relabeling(result))
    ordered_slots = sortperm(1:nloops; by=slot -> rank[pair_vertices[slot]])
    loop_permutation = zeros(Int, nloops)
    loop_signs = ones(Int, nloops)
    for (new_slot, old_slot) in enumerate(ordered_slots)
        loop_permutation[old_slot] = new_slot
        incidences = loop_incidences[old_slot]
        isempty(incidences) && continue
        canonical_first = first(sort!(copy(incidences); by=record -> rank[first(record)]))
        loop_signs[old_slot] = last(canonical_first) > 0 ? 1 : -1
    end
    return KC.loop_permutation_transform(
        basis, external, loop_permutation, loop_signs, zeros(Int, nloops)
    )
end

function gc_quotient_loop_momenta(
    expression::KC.OccupationReducedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:KC.Statistics,O,G,Ctx<:KC.AbstractWignerContext}
    D = promote_type(C, KC.ComplexRationals, Rational{Int})
    out = Dict{KC.CollisionKernelSector{S},KC.OccupationPolynomial{D,S}}()

    for (sector, polynomial) in KC.occupation_reduced_terms(expression)
        for (occupation_monomial, occupation_coefficient) in polynomial
            for (kinematic_monomial, kinematic_coefficient) in KC.kinematic_factor(sector)
                atom_sector = KC._kinematic_atom_sector(sector, kinematic_monomial)
                transform = gc_canonical_loop_transform(atom_sector, occupation_monomial)
                transformed_sector, support_factor = KC._transform_kernel_sector(
                    atom_sector, transform
                )
                transformed_monomial = KC.transform_loop_momenta(
                    occupation_monomial, transform
                )
                transformed_coefficient =
                    convert(D, occupation_coefficient) *
                    convert(D, kinematic_coefficient) *
                    convert(D, support_factor)
                contribution = Pair{KC.OccupationMonomial{S},D}[transformed_monomial => transformed_coefficient]
                KC._push_kernel_polynomial!(
                    out, transformed_sector, KC.OccupationPolynomial{D,S}(contribution)
                )
            end
        end
    end

    return KC.LoopQuotientedExpression{D,S,O,G,Ctx}(
        out,
        KC.target_family(expression),
        KC.parameters(expression),
        KC.wigner_context(expression),
    )
end

function assert_quotient_matches(expression)
    nauty = KC.quotient_loop_momenta(expression)
    gc = gc_quotient_loop_momenta(expression)
    @test KC.loop_quotient_terms(gc) == KC.loop_quotient_terms(nauty)
    @test KC.collision_kernel(gc).terms == KC.collision_kernel(nauty).terms
    return nothing
end

function benchmark_pair(label, nauty, gc)
    nauty_trial = @benchmark $nauty() samples = 25 evals = 1
    gc_trial = @benchmark $gc() samples = 25 evals = 1
    nauty_estimate = median(nauty_trial)
    gc_estimate = median(gc_trial)
    println(
        label,
        ": Nauty ",
        round(nauty_estimate.time / 1.0e3; digits=2),
        " μs / ",
        nauty_estimate.memory,
        " B / ",
        nauty_estimate.allocs,
        " allocs; GC ",
        round(gc_estimate.time / 1.0e3; digits=2),
        " μs / ",
        gc_estimate.memory,
        " B / ",
        gc_estimate.allocs,
        " allocs",
    )
    return nothing
end

@qfields gc_adapter_ϕ::Boson
@qfields gc_adapter_χ::Boson

c, q = gc_adapter_ϕ[Classical], gc_adapter_ϕ[Quantum]
χc, χq = gc_adapter_χ[Classical], gc_adapter_χ[Quantum]

function as_contractions(vs)
    return KC.Contraction{Boson}[KC.Contraction(item) for item in vs]
end

@testset "GraphCombinations exact KC oracle adapter ($GC_SHA)" begin
    @testset "propagator canonicalization and automorphism order" begin
        ring = as_contractions([
            (c(Out()), bar(q)(Bulk(1))),
            (c(Bulk(1)), bar(q)(Bulk(2))),
            (c(Bulk(2)), bar(q)(Bulk(3))),
            (c(Bulk(3)), bar(q)(Bulk(4))),
            (c(Bulk(4)), bar(q)(Bulk(1))),
            (c(Bulk(4)), bar(q)(In())),
        ])
        self_loop = as_contractions([
            (c(Out()), bar(q)(Bulk(1))),
            (c(Bulk(1)), bar(q)(Bulk(1))),
            (c(Bulk(1)), bar(q)(In())),
        ])
        repeated = as_contractions([
            (c(Out()), bar(q)(Bulk(1))),
            (c(Bulk(1)), bar(q)(Bulk(2))),
            (c(Bulk(1)), bar(q)(Bulk(2))),
            (c(Bulk(2)), bar(q)(In())),
        ])
        colored = as_contractions([
            (c(Out()), bar(q)(Bulk(1))),
            (χc(Bulk(1)), bar(χq)(Bulk(2))),
            (c(Bulk(2)), bar(q)(In())),
        ])
        symmetric = as_contractions([
            (c(Bulk(1)), bar(q)(Bulk(2))),
            (c(Bulk(2)), bar(q)(Bulk(3))),
            (c(Bulk(3)), bar(q)(Bulk(1))),
        ])

        for fixture in (ring, self_loop, repeated, colored, symmetric)
            gc_fixture = gc_canonicalize(fixture)
            nauty_fixture = KC.canonicalize(fixture)
            @test gc_canonicalize(nauty_fixture) == gc_fixture
            @test KC.canonicalize(gc_fixture) == nauty_fixture

            graph_positions = KC.canonicalization_positions(fixture)
            gc_topology = gc_topology_result(fixture, graph_positions)
            _, _, _, nauty_automorphisms = KC.canonicalization_permutations(
                fixture, graph_positions
            )
            @test GC.canonical_automorphism_order(gc_topology) == nauty_automorphisms.n
        end
        @test GC.canonical_automorphism_order(gc_topology_result(symmetric)) == 3
    end

    @testset "historical topology 1 / 3 / 11 / 59" begin
        elastic = -(
            1 // 2 * (c^2 + q^2) * bar(c) * bar(q) + 1 // 2 * c * q * (bar(c)^2 + bar(q)^2)
        )
        L = InteractionLagrangian(elastic)
        cases = ((1, 3, 1), (2, 5, 3), (3, 7, 11), (4, 9, 59))
        for (order, edge_count, expected) in cases
            G = DressedPropagator(L, Val(order), Val(edge_count))
            component = KC.topologies(G.keldysh)
            @test length(keys(component)) == expected

            nauty_to_gc = Dict{Any,Any}()
            gc_to_nauty = Dict{Any,Any}()
            for (key, diagrams) in component
                nauty_key = Tuple(key)
                for diagram in diagrams
                    contractions = KC.Contraction{Boson}[
                        (edge.out, edge.in) for edge in KC.contractions(diagram)
                    ]
                    gc_key = Tuple(gc_legacy_topology(contractions, Val(length(key))))
                    @test get!(nauty_to_gc, nauty_key, gc_key) == gc_key
                    @test get!(gc_to_nauty, gc_key, nauty_key) == nauty_key
                end
            end
            @test length(nauty_to_gc) == expected
            @test length(gc_to_nauty) == expected
        end
    end

    include(joinpath(@__DIR__, "..", "benchmarks", "collision_reduction.jl"))
    include(joinpath(@__DIR__, "..", "benchmarks", "fermionic_pwave_loss.jl"))

    @testset "loop-momentum quotient fixtures" begin
        for nloops in (2, 4)
            expression = benchmark_loop_quotient_fixture(nloops)
            assert_quotient_matches(expression)
            for (sector, polynomial) in KC.occupation_reduced_terms(expression)
                for (monomial, _) in polynomial
                    for (kinematic, _) in KC.kinematic_factor(sector)
                        atom_sector = KC._kinematic_atom_sector(sector, kinematic)
                        gc_transform = gc_canonical_loop_transform(atom_sector, monomial)
                        nauty_transform = KC._canonical_loop_transform(
                            atom_sector, monomial
                        )
                        @test KC.loop_transform_matrix(gc_transform) ==
                            KC.loop_transform_matrix(nauty_transform)
                    end
                end
            end
        end
    end

    @testset "production collision cases" begin
        bosonic_collision = benchmark_collision_fixture()
        bosonic_reduced = KC.reduce_frequency_collision(bosonic_collision)
        bosonic_occupation = KC.occupation_reduced_expression(bosonic_reduced)
        assert_quotient_matches(bosonic_occupation)

        _, _, _, _, _, _, _, _, _, fermionic_occupation, fermionic_quotient = benchmark_fermionic_pwave_fixtures()
        gc_fermionic_quotient = gc_quotient_loop_momenta(fermionic_occupation)
        @test KC.loop_quotient_terms(gc_fermionic_quotient) ==
            KC.loop_quotient_terms(fermionic_quotient)
        @test KC.collision_kernel(gc_fermionic_quotient).terms ==
            KC.collision_kernel(fermionic_quotient).terms
    end

    @testset "backend measurements" begin
        ring = as_contractions([
            (c(Out()), bar(q)(Bulk(1))),
            (c(Bulk(1)), bar(q)(Bulk(2))),
            (c(Bulk(2)), bar(q)(Bulk(3))),
            (c(Bulk(3)), bar(q)(Bulk(4))),
            (c(Bulk(4)), bar(q)(Bulk(1))),
            (c(Bulk(4)), bar(q)(In())),
        ])
        benchmark_pair(
            "propagator ring", () -> KC.canonicalize(ring), () -> gc_canonicalize(ring)
        )

        occupation2 = benchmark_loop_quotient_fixture(2)
        occupation4 = benchmark_loop_quotient_fixture(4)
        benchmark_pair(
            "loop quotient 2-loop",
            () -> KC.quotient_loop_momenta(occupation2),
            () -> gc_quotient_loop_momenta(occupation2),
        )
        benchmark_pair(
            "loop quotient 4-loop",
            () -> KC.quotient_loop_momenta(occupation4),
            () -> gc_quotient_loop_momenta(occupation4),
        )
    end
end
