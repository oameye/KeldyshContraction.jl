using Test
using Combinatorics

import GraphCombinations as GC
import KeldyshContraction as KC
using KeldyshContraction: @qfields, Boson

function _gc_direct_graph(vs, graph_positions)
    edges = Pair{Int,Int}[]
    seen = Set{Tuple{Int,Int}}()
    simple = true
    for item in vs
        source_position, target_position = KC.positions(item)
        source = KC.position_vertex(graph_positions, source_position)
        target = KC.position_vertex(graph_positions, target_position)
        key = (source, target)
        if key in seen
            simple = false
        else
            push!(seen, key)
            push!(edges, source => target)
        end
    end
    return GC.DirectedGCGraph(edges, length(graph_positions)), KC.position_labels(graph_positions), simple
end

function _gc_colored_graph(vs, graph_positions)
    colors = KC.propagator_colors(vs)
    npositions = length(graph_positions)
    labels = Vector{Int}(undef, npositions + length(vs))
    copyto!(labels, 1, KC.position_labels(graph_positions), 1, npositions)

    edges = Pair{Int,Int}[]
    sizehint!(edges, 2 * length(vs))
    for (i, item) in enumerate(vs)
        color_index = searchsortedfirst(colors, KC.propagator_color(item))
        edge_vertex = npositions + i
        labels[edge_vertex] = 3 + color_index
        source_position, target_position = KC.positions(item)
        push!(
            edges,
            KC.position_vertex(graph_positions, source_position) => edge_vertex,
        )
        push!(
            edges,
            edge_vertex => KC.position_vertex(graph_positions, target_position),
        )
    end
    return GC.DirectedGCGraph(edges, length(labels)), labels
end

function gc_physical_permutation(vs, graph_positions)
    direct_graph, direct_colors, simple = _gc_direct_graph(vs, graph_positions)
    capacity = length(graph_positions) + length(vs)
    workspace = GC.DirectedCanonicalizationWorkspace(capacity)
    buffer = GC.DirectedCanonicalizationBuffer(capacity)
    GC.canonicalize_directed!(buffer, workspace, direct_graph, direct_colors)

    use_direct = isone(GC.canonical_automorphism_order(buffer)) ||
                 (simple && KC.uniform_coloring(vs))
    if !use_direct
        colored_graph, colored_colors = _gc_colored_graph(vs, graph_positions)
        GC.canonicalize_directed!(buffer, workspace, colored_graph, colored_colors)
    end
    return [GC.canonical_rank(buffer, old_vertex) for old_vertex in 1:buffer.num_vertices]
end

function gc_physical_canonicalize(vs::Vector{T}) where {T}
    isempty(vs) && return copy(vs)
    graph_positions = KC.canonicalization_positions(vs)
    permutation = gc_physical_permutation(vs, graph_positions)
    mapping = KC.make_permutation_dict(permutation, graph_positions, vs)
    return T[KC.relabel_bulk_positions(item, mapping) for item in vs]
end

function bulk_count(vs)
    positions = KC.canonicalization_positions(vs)
    return count(KC.is_bulk, positions)
end

function permute_bulk_labels(vs::Vector{T}, permutation) where {T}
    mapping = Dict{KC.Position,KC.Position}(
        KC.Bulk(i) => KC.Bulk(permutation[i]) for i in eachindex(permutation)
    )
    return T[KC.relabel_bulk_positions(item, mapping) for item in vs]
end

function add_orbit!(cases, label, vs; exhaustive=true)
    push!(cases, label => vs)
    nbulk = bulk_count(vs)
    if exhaustive && nbulk > 1
        for (index, permutation) in enumerate(permutations(1:nbulk))
            permutation == collect(1:nbulk) && continue
            push!(cases, "$label/relabel-$index" => permute_bulk_labels(vs, permutation))
        end
    end
    return nothing
end

@qfields ϕ::Boson
c, q = ϕ[KC.Classical], ϕ[KC.Quantum]
@qfields χ::Boson
χc, χq = χ[KC.Classical], χ[KC.Quantum]

cases = Pair{String,Any}[]

add_orbit!(
    cases,
    "linear-2-bulk",
    [
        KC.Contraction(c(KC.Out()), KC.bar(q)(KC.Bulk(1))),
        KC.Contraction(c(KC.Bulk(1)), KC.bar(q)(KC.Bulk(2))),
        KC.Contraction(c(KC.Bulk(2)), KC.bar(q)(KC.In())),
    ],
)

add_orbit!(
    cases,
    "colored-linear",
    [
        KC.Contraction(c(KC.Out()), KC.bar(q)(KC.Bulk(1))),
        KC.Contraction(χc(KC.Bulk(1)), KC.bar(χq)(KC.Bulk(2))),
        KC.Contraction(c(KC.Bulk(2)), KC.bar(q)(KC.In())),
    ],
)

add_orbit!(
    cases,
    "self-loop",
    [
        KC.Contraction(c(KC.Out()), KC.bar(q)(KC.Bulk(1))),
        KC.Contraction(c(KC.Bulk(1)), KC.bar(q)(KC.Bulk(1))),
        KC.Contraction(c(KC.Bulk(1)), KC.bar(q)(KC.In())),
    ],
)

add_orbit!(
    cases,
    "parallel-edge",
    [
        KC.Contraction(c(KC.Out()), KC.bar(q)(KC.Bulk(1))),
        KC.Contraction(c(KC.Bulk(1)), KC.bar(q)(KC.Bulk(2))),
        KC.Contraction(c(KC.Bulk(1)), KC.bar(q)(KC.Bulk(2))),
        KC.Contraction(c(KC.Bulk(2)), KC.bar(q)(KC.In())),
    ],
)

add_orbit!(
    cases,
    "three-bulk",
    [
        KC.Contraction(c(KC.Out()), KC.bar(q)(KC.Bulk(1))),
        KC.Contraction(c(KC.Bulk(1)), KC.bar(q)(KC.Bulk(2))),
        KC.Contraction(c(KC.Bulk(1)), KC.bar(q)(KC.Bulk(3))),
        KC.Contraction(c(KC.Bulk(3)), KC.bar(q)(KC.Bulk(1))),
        KC.Contraction(c(KC.Bulk(2)), KC.bar(q)(KC.Bulk(3))),
        KC.Contraction(c(KC.Bulk(3)), KC.bar(q)(KC.Bulk(2))),
        KC.Contraction(c(KC.Bulk(2)), KC.bar(q)(KC.In())),
    ],
)

add_orbit!(
    cases,
    "ring-4-bulk",
    [
        KC.Contraction(c(KC.Out()), KC.bar(q)(KC.Bulk(1))),
        KC.Contraction(c(KC.Bulk(1)), KC.bar(q)(KC.Bulk(2))),
        KC.Contraction(c(KC.Bulk(2)), KC.bar(q)(KC.Bulk(3))),
        KC.Contraction(c(KC.Bulk(3)), KC.bar(q)(KC.Bulk(4))),
        KC.Contraction(c(KC.Bulk(4)), KC.bar(q)(KC.Bulk(1))),
        KC.Contraction(c(KC.Bulk(4)), KC.bar(q)(KC.In())),
    ],
)

add_orbit!(
    cases,
    "third-order-two-body",
    [
        KC.Contraction(c(KC.Out()), KC.bar(q)(KC.Bulk(1))),
        KC.Contraction(c(KC.Bulk(1)), KC.bar(q)(KC.Bulk(2))),
        KC.Contraction(c(KC.Bulk(2)), KC.bar(q)(KC.Bulk(2))),
        KC.Contraction(c(KC.Bulk(3)), KC.bar(q)(KC.Bulk(1))),
        KC.Contraction(c(KC.Bulk(2)), KC.bar(q)(KC.Bulk(3))),
        KC.Contraction(c(KC.Bulk(3)), KC.bar(q)(KC.Bulk(3))),
        KC.Contraction(c(KC.Bulk(1)), KC.bar(q)(KC.In())),
    ],
)

println("GC/Nauty physical canonicalization oracle: ", length(cases), " cases")
permutation_mismatches = String[]
physical_mismatches = String[]

@testset "GC physical canonicalization oracle" begin
    for (label, vs) in cases
        graph_positions = KC.canonicalization_positions(vs)
        gc_perm = gc_physical_permutation(vs, graph_positions)
        nauty_perm = KC.canonicalization_permutation(vs, graph_positions)
        if gc_perm != nauty_perm
            push!(permutation_mismatches, label)
        end

        gc_canonical = gc_physical_canonicalize(vs)
        nauty_canonical = KC.canonicalize(vs)
        if gc_canonical != nauty_canonical
            push!(physical_mismatches, label)
        end
        @test gc_canonical == nauty_canonical
    end
end

println("canonical-permutation mismatches: ", length(permutation_mismatches))
for label in permutation_mismatches
    println("  ", label)
end
println("physical canonical-output mismatches: ", length(physical_mismatches))
for label in physical_mismatches
    println("  ", label)
end
