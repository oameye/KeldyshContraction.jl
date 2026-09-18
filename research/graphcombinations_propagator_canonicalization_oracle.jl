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
        push!(edges, KC.position_vertex(graph_positions, source_position) => edge_vertex)
        push!(edges, edge_vertex => KC.position_vertex(graph_positions, target_position))
    end
    return GC.DirectedGCGraph(edges, length(labels)), labels
end

function gc_physical_witnesses(vs, graph_positions)
    direct_graph, direct_colors, simple = _gc_direct_graph(vs, graph_positions)
    capacity = length(graph_positions) + length(vs)
    workspace = GC.DirectedCanonicalizationWorkspace(capacity)
    buffer = GC.DirectedCanonicalizationBuffer(capacity)
    GC.canonicalize_directed!(buffer, workspace, direct_graph, direct_colors)

    direct_automorphisms = GC.canonical_automorphism_order(buffer)
    use_direct = isone(direct_automorphisms) || (simple && KC.uniform_coloring(vs))
    if !use_direct
        colored_graph, colored_colors = _gc_colored_graph(vs, graph_positions)
        GC.canonicalize_directed!(buffer, workspace, colored_graph, colored_colors)
    end

    n = buffer.num_vertices
    canonical_to_old = [GC.original_vertex(buffer, rank) for rank in 1:n]
    old_to_canonical = [GC.canonical_rank(buffer, old_vertex) for old_vertex in 1:n]
    return (;
        canonical_to_old,
        old_to_canonical,
        use_direct,
        simple,
        direct_automorphisms,
        final_automorphisms=GC.canonical_automorphism_order(buffer),
    )
end

function canonicalize_with_permutation(vs::Vector{T}, graph_positions, permutation) where {T}
    mapping = KC.make_permutation_dict(permutation, graph_positions, vs)
    canonical = T[KC.relabel_bulk_positions(item, mapping) for item in vs]
    return canonical, mapping
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

family(label) = first(split(label, "/"))

function compact_map(mapping)
    pairs = collect(mapping)
    sort!(pairs; by=p -> string(first(p)))
    return join(("$(first(p))=>$(last(p))" for p in pairs), ", ")
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

println("GC/Nauty physical canonicalization diagnostic: ", length(cases), " cases")

records = NamedTuple[]
seed_gc_outputs = Dict{String,Any}()
seed_nauty_outputs = Dict{String,Any}()

@testset "GC witness consistency" begin
    for (label, vs) in cases
        graph_positions = KC.canonicalization_positions(vs)
        witnesses = gc_physical_witnesses(vs, graph_positions)
        n = length(witnesses.canonical_to_old)
        @test sort(witnesses.canonical_to_old) == collect(1:n)
        @test sort(witnesses.old_to_canonical) == collect(1:n)
        @test all(
            witnesses.canonical_to_old[witnesses.old_to_canonical[old_vertex]] == old_vertex for
            old_vertex in 1:n
        )

        nauty_perm = KC.canonicalization_permutation(vs, graph_positions)
        gc_original, gc_original_map =
            canonicalize_with_permutation(vs, graph_positions, witnesses.canonical_to_old)
        gc_inverse, gc_inverse_map =
            canonicalize_with_permutation(vs, graph_positions, witnesses.old_to_canonical)
        nauty_output, nauty_map = canonicalize_with_permutation(vs, graph_positions, nauty_perm)
        @test nauty_output == KC.canonicalize(vs)

        fam = family(label)
        if label == fam
            seed_gc_outputs[fam] = gc_original
            seed_nauty_outputs[fam] = nauty_output
        end

        push!(records, (;
            label,
            family=fam,
            graph_positions,
            witnesses,
            nauty_perm,
            gc_original_map,
            gc_inverse_map,
            nauty_map,
            original_perm_equal=witnesses.canonical_to_old == nauty_perm,
            inverse_perm_equal=witnesses.old_to_canonical == nauty_perm,
            original_map_equal=gc_original_map == nauty_map,
            inverse_map_equal=gc_inverse_map == nauty_map,
            original_output_equal=gc_original == nauty_output,
            inverse_output_equal=gc_inverse == nauty_output,
            gc_original,
            gc_inverse,
            nauty_output,
        ))
    end
end

println()
println("aggregate comparison")
for field in (
    :original_perm_equal,
    :inverse_perm_equal,
    :original_map_equal,
    :inverse_map_equal,
    :original_output_equal,
    :inverse_output_equal,
)
    count_equal = count(r -> getproperty(r, field), records)
    println("  ", field, ": ", count_equal, "/", length(records))
end

println()
println("by family")
for fam in unique(r.family for r in records)
    rs = filter(r -> r.family == fam, records)
    direct = count(r -> r.witnesses.use_direct, rs)
    original = count(r -> r.original_output_equal, rs)
    inverse = count(r -> r.inverse_output_equal, rs)
    original_maps = count(r -> r.original_map_equal, rs)
    inverse_maps = count(r -> r.inverse_map_equal, rs)
    aut_orders = sort!(unique(r.witnesses.final_automorphisms for r in rs))
    println(
        "  ", fam,
        ": cases=", length(rs),
        " direct=", direct,
        " original-output=", original,
        " inverse-output=", inverse,
        " original-map=", original_maps,
        " inverse-map=", inverse_maps,
        " final-aut=", aut_orders,
    )
end

println()
println("relabeling invariance against each family's seed")
for fam in unique(r.family for r in records)
    rs = filter(r -> r.family == fam, records)
    gc_seed = seed_gc_outputs[fam]
    nauty_seed = seed_nauty_outputs[fam]
    gc_invariant = count(r -> r.gc_original == gc_seed, rs)
    nauty_invariant = count(r -> r.nauty_output == nauty_seed, rs)
    println(
        "  ", fam,
        ": GC-original=", gc_invariant, "/", length(rs),
        " Nauty=", nauty_invariant, "/", length(rs),
    )
end

mismatches = filter(r -> !r.original_output_equal, records)
println()
println("first detailed original-witness mismatches: ", min(length(mismatches), 12))
for r in Iterators.take(mismatches, 12)
    w = r.witnesses
    println("CASE ", r.label)
    println(
        "  path=", w.use_direct ? "direct" : "colored",
        " simple=", w.simple,
        " direct-aut=", w.direct_automorphisms,
        " final-aut=", w.final_automorphisms,
    )
    println("  positions=", r.graph_positions)
    println("  nauty canonical->old=", r.nauty_perm)
    println("  GC canonical->old=", w.canonical_to_old)
    println("  GC old->canonical=", w.old_to_canonical)
    println("  nauty map={", compact_map(r.nauty_map), "}")
    println("  GC original map={", compact_map(r.gc_original_map), "}")
    println("  GC inverse map={", compact_map(r.gc_inverse_map), "}")
    println(
        "  map equality: original=", r.original_map_equal,
        " inverse=", r.inverse_map_equal,
        "; output equality: original=", r.original_output_equal,
        " inverse=", r.inverse_output_equal,
    )
end
