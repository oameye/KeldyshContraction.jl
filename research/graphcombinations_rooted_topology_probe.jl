using KeldyshContraction
using Test

using Combinatorics: Combinatorics
import GraphCombinations as GC
import KeldyshContraction as KC

const GC_SHA = "e0e5d7dc4e14f4fc356ffeee3809fdb5b47fcf12"

function gc_result(edges::Vector{Tuple{Int,Int}}, colors::Vector{Int})
    graph = GC.DirectedGCGraph(
        Pair{Int,Int}[source => target for (source, target) in edges], length(colors)
    )
    return graph, GC.canonicalize_directed(graph, colors)
end

function rooted_graph_key(
    vs;
    directed::Bool,
    preserve_multiplicity::Bool,
    include_loops::Bool,
    colored_roots::Bool,
)
    graph_positions = KC.canonicalization_positions(vs)
    edges = Tuple{Int,Int}[]
    for item in vs
        out, in = KC.positions(item)
        source = KC.position_vertex(graph_positions, out)
        target = KC.position_vertex(graph_positions, in)
        !include_loops && source == target && continue
        push!(edges, (source, target))
        !directed && source != target && push!(edges, (target, source))
    end
    preserve_multiplicity || unique!(edges)

    colors = colored_roots ? KC.position_labels(graph_positions) : ones(Int, length(graph_positions))
    graph, result = gc_result(edges, colors)
    old_to_canonical = GC.vertex_mapping(GC.canonical_relabeling(result))
    canonical_colors = similar(colors)
    for old_vertex in eachindex(colors)
        canonical_colors[old_to_canonical[old_vertex]] = colors[old_vertex]
    end

    buffer = GC.DirectedCanonicalizationBuffer(graph.num_vertices)
    workspace = GC.DirectedCanonicalizationWorkspace(graph.num_vertices)
    GC.canonicalize_directed!(buffer, workspace, graph, colors)
    @test Tuple(buffer.canonical_multiplicities) ==
          Tuple(GC.canonical_graph(result).multiplicities)
    @test buffer.old_to_canonical == old_to_canonical

    return (
        Tuple(canonical_colors), Tuple(GC.canonical_graph(result).multiplicities)
    )
end

function rooted_primary_signature(vs, mapping, nbulk::Int)
    nvertices = nbulk + 2
    adjacency = zeros(UInt8, nvertices * nvertices)
    for item in vs
        out, in = KC.positions(item)
        source = if KC.is_out(out)
            1
        elseif KC.is_in(out)
            nvertices
        else
            1 + mapping[out]
        end
        target = if KC.is_out(in)
            1
        elseif KC.is_in(in)
            nvertices
        else
            1 + mapping[in]
        end
        adjacency[(source - 1) * nvertices + target] = 1
    end
    return Tuple(adjacency)
end

function two_stage_rooted_topology(vs, ::Val{E2}) where {E2}
    graph_positions = KC.canonicalization_positions(vs)
    bulk_positions = KC.Position[p for p in graph_positions if KC.is_bulk(p)]
    nbulk = length(bulk_positions)
    best_primary = nothing
    best_topology = nothing

    for bulk_order in Combinatorics.permutations(bulk_positions)
        rank = Dict(position => i for (i, position) in enumerate(bulk_order))
        primary = rooted_primary_signature(vs, rank, nbulk)
        position_mapping = Dict(position => KC.Bulk(i) for (i, position) in enumerate(bulk_order))
        topology_edges = Tuple{Int8,Int8}[
            KC.integer_positions(KC.relabel_bulk_positions(item, position_mapping)) for item in vs
        ]
        topology = Tuple(KC.bulk_multiplicity(topology_edges, Val(E2)))

        if best_primary === nothing || primary < best_primary
            best_primary = primary
            best_topology = topology
        elseif primary == best_primary && topology < best_topology
            best_topology = topology
        end
    end
    return best_topology
end

const ROOTED_VARIANTS = (
    historical_image=(
        directed=true,
        preserve_multiplicity=false,
        include_loops=true,
        colored_roots=false,
    ),
    rooted_directed_simple=(
        directed=true,
        preserve_multiplicity=false,
        include_loops=true,
        colored_roots=true,
    ),
    rooted_directed_multi=(
        directed=true,
        preserve_multiplicity=true,
        include_loops=true,
        colored_roots=true,
    ),
    rooted_undirected_multi=(
        directed=false,
        preserve_multiplicity=true,
        include_loops=true,
        colored_roots=true,
    ),
    rooted_undirected_multi_no_loops=(
        directed=false,
        preserve_multiplicity=true,
        include_loops=false,
        colored_roots=true,
    ),
)

@qfields rooted_probe_ϕ::Boson
c, q = rooted_probe_ϕ[Classical], rooted_probe_ϕ[Quantum]

println("GC rooted topology probe ($GC_SHA)")
flush(stdout)

@testset "rooted topology semantic graph diagnostics" begin
    elastic = -(
        1 // 2 * (c^2 + q^2) * bar(c) * bar(q) + 1 // 2 * c * q * (bar(c)^2 + bar(q)^2)
    )
    L = InteractionLagrangian(elastic)

    # Order three is the first historical mismatch. Use it to identify the minimal
    # backend-neutral rooted graph before paying for the order-four certification.
    order, edge_count, expected = (3, 7, 11)
    G = DressedPropagator(L, Val(order), Val(edge_count))
    component = KC.topologies(G.keldysh)
    @test length(keys(component)) == expected

    for (variant, options) in pairs(ROOTED_VARIANTS)
        semantic_to_nauty = Dict{Any,Any}()
        within_class_splits = 0
        cross_class_collisions = 0

        for (key, diagrams) in component
            nauty_key = Tuple(key)
            semantic_keys = Set{Any}()
            for diagram in diagrams
                contractions = KC.Contraction{Boson}[
                    (edge.out, edge.in) for edge in KC.contractions(diagram)
                ]
                push!(semantic_keys, rooted_graph_key(contractions; options...))
            end

            within_class_splits += length(semantic_keys) != 1
            for semantic_key in semantic_keys
                previous = get(semantic_to_nauty, semantic_key, nothing)
                cross_class_collisions += previous !== nothing && previous != nauty_key
                semantic_to_nauty[semantic_key] = nauty_key
            end
        end

        println(
            variant,
            ": keys=",
            length(semantic_to_nauty),
            ", within-class splits=",
            within_class_splits,
            ", cross-class collisions=",
            cross_class_collisions,
        )
        flush(stdout)
    end

    topology_to_nauty = Dict{Any,Any}()
    nauty_to_topology = Dict{Any,Any}()
    for (key, diagrams) in component
        nauty_key = Tuple(key)
        topology_keys = Set{Any}()
        for diagram in diagrams
            contractions = KC.Contraction{Boson}[
                (edge.out, edge.in) for edge in KC.contractions(diagram)
            ]
            push!(topology_keys, two_stage_rooted_topology(contractions, Val(length(key))))
        end
        @test length(topology_keys) == 1
        topology_key = only(topology_keys)
        @test get!(nauty_to_topology, nauty_key, topology_key) == topology_key
        previous = get(topology_to_nauty, topology_key, nothing)
        @test previous === nothing || previous == nauty_key
        topology_to_nauty[topology_key] = nauty_key
    end
    println("two-stage rooted topology: ", length(topology_to_nauty), " classes")
    @test length(topology_to_nauty) == expected
    @test length(nauty_to_topology) == expected
end
