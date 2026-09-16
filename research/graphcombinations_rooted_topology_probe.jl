using KeldyshContraction
using Test

import GraphCombinations as GC
import KeldyshContraction as KC

const GC_SHA = "e0e5d7dc4e14f4fc356ffeee3809fdb5b47fcf12"

function gc_result(edges::Vector{Tuple{Int,Int}}, num_vertices::Int)
    graph = GC.DirectedGCGraph(
        Pair{Int,Int}[source => target for (source, target) in edges], num_vertices
    )
    return graph, GC.canonicalize_directed(graph)
end

function legacy_graph_data(vs)
    position_pairs = Tuple{Int8,Int8}[KC.integer_positions(item) for item in vs]
    flattened = collect(Iterators.flatten(position_pairs))
    max_label = length(unique(flattened))
    has_out = typemin(Int8) in flattened

    edges = Tuple{Int,Int}[]
    sizehint!(edges, length(position_pairs))
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
    return edges, max_label, has_out
end

function gc_legacy_witness_topology(vs, ::Val{E2}) where {E2}
    isempty(vs) && return Tuple(KC.bulk_multiplicity(Tuple{Int8,Int8}[], Val(E2)))

    edges, max_label, has_out = legacy_graph_data(vs)
    _, result = gc_result(edges, max_label)
    old_to_canonical = GC.vertex_mapping(GC.canonical_relabeling(result))
    canonical_to_old = sortperm(old_to_canonical)
    mapping = KC.legacy_topology_permutation_dict(canonical_to_old, max_label, has_out)
    topology_edges = Tuple{Int8,Int8}[
        KC.integer_positions(KC.relabel_bulk_positions(item, mapping)) for item in vs
    ]
    return Tuple(KC.bulk_multiplicity(topology_edges, Val(E2)))
end

function gc_semantic_topology(vs)
    position_pairs = Tuple{Int8,Int8}[KC.integer_positions(item) for item in vs]
    bulk_labels = sort!(
        unique(
            Int(label) for pair in position_pairs for label in pair if
            label != typemin(Int8) && label != typemax(Int8)
        ),
    )
    label_to_vertex = Dict(label => i for (i, label) in enumerate(bulk_labels))

    edges = Tuple{Int,Int}[]
    for pair in position_pairs
        a, b = pair
        if a == typemin(Int8) ||
            a == typemax(Int8) ||
            b == typemin(Int8) ||
            b == typemax(Int8) ||
            a == b
            continue
        end
        u = label_to_vertex[Int(a)]
        v = label_to_vertex[Int(b)]
        push!(edges, (u, v))
        push!(edges, (v, u))
    end

    graph, result = gc_result(edges, length(bulk_labels))
    canonical = Tuple(GC.canonical_graph(result).multiplicities)

    buffer = GC.DirectedCanonicalizationBuffer(graph.num_vertices)
    workspace = GC.DirectedCanonicalizationWorkspace(graph.num_vertices)
    GC.canonicalize_directed!(buffer, workspace, graph)
    @test Tuple(buffer.canonical_multiplicities) == canonical
    @test buffer.old_to_canonical ==
          GC.vertex_mapping(GC.canonical_relabeling(result))

    return canonical
end

@qfields rooted_probe_ϕ::Boson
c, q = rooted_probe_ϕ[Classical], rooted_probe_ϕ[Quantum]

println("GC rooted topology probe ($GC_SHA)")
flush(stdout)

@testset "semantic GC topology partition versus historical Nauty partition" begin
    elastic = -(
        1 // 2 * (c^2 + q^2) * bar(c) * bar(q) + 1 // 2 * c * q * (bar(c)^2 + bar(q)^2)
    )
    L = InteractionLagrangian(elastic)
    cases = ((1, 3, 1), (2, 5, 3), (3, 7, 11), (4, 9, 59))

    for (order, edge_count, expected) in cases
        G = DressedPropagator(L, Val(order), Val(edge_count))
        component = KC.topologies(G.keldysh)
        @test length(keys(component)) == expected

        semantic_to_nauty = Dict{Any,Any}()
        nauty_to_semantic = Dict{Any,Any}()
        first_legacy_split = nothing
        first_legacy_collapse = nothing
        legacy_to_nauty = Dict{Any,Any}()

        for (key, diagrams) in component
            nauty_key = Tuple(key)
            semantic_keys = Set{Any}()
            legacy_keys = Set{Any}()

            for diagram in diagrams
                contractions = KC.Contraction{Boson}[
                    (edge.out, edge.in) for edge in KC.contractions(diagram)
                ]
                push!(semantic_keys, gc_semantic_topology(contractions))
                push!(legacy_keys, gc_legacy_witness_topology(contractions, Val(length(key))))
            end

            @test length(semantic_keys) == 1
            semantic_key = only(semantic_keys)
            @test get!(nauty_to_semantic, nauty_key, semantic_key) == semantic_key
            previous = get(semantic_to_nauty, semantic_key, nothing)
            @test previous === nothing || previous == nauty_key
            semantic_to_nauty[semantic_key] = nauty_key

            if length(legacy_keys) != 1 && first_legacy_split === nothing
                first_legacy_split = (nauty_key=nauty_key, legacy_keys=collect(legacy_keys))
            elseif length(legacy_keys) == 1
                legacy_key = only(legacy_keys)
                previous_legacy = get(legacy_to_nauty, legacy_key, nothing)
                if previous_legacy !== nothing &&
                    previous_legacy != nauty_key &&
                    first_legacy_collapse === nothing
                    first_legacy_collapse = (
                        first_nauty_key=previous_legacy,
                        second_nauty_key=nauty_key,
                        legacy_key=legacy_key,
                    )
                end
                legacy_to_nauty[legacy_key] = nauty_key
            end
        end

        if first_legacy_split !== nothing
            println("legacy witness split at order $order: $first_legacy_split")
        end
        if first_legacy_collapse !== nothing
            println("legacy witness collapse at order $order: $first_legacy_collapse")
        end
        flush(stdout)

        @test length(nauty_to_semantic) == expected
        @test length(semantic_to_nauty) == expected
        println("semantic topology order $order -> $expected classes")
        flush(stdout)
    end
end
