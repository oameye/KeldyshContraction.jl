using KeldyshContraction, Test
using Combinatorics: permutations
import KeldyshContraction as KC
const GraphComb = KC.GC

function signed_pairing_weights(pairings)
    weights = Dict{Any,Int}()
    for (pairing, topology, multiplicity) in pairings
        key = (pairing.contractions, topology)
        weight = Int(pairing.sign) * multiplicity
        weights[key] = get(weights, key, 0) + weight
    end
    filter!(pair -> !iszero(last(pair)), weights)
    return weights
end

function port_cells(indices, ncolors)
    vertices = Vector{Int}(undef, length(indices))
    colors = Vector{Int}(undef, length(indices))
    for i in eachindex(indices)
        index = indices[i]
        vertices[i] = fld(index - 1, ncolors) + 1
        colors[i] = mod(index - 1, ncolors) + 1
    end
    return vertices, colors
end

function port_counts(vertices, colors, nvertices, ncolors)
    counts = zeros(Int, nvertices, ncolors)
    for i in eachindex(vertices)
        counts[vertices[i], colors[i]] += 1
    end
    return counts
end

function canonical_port_edge_key(problem, edges, nvertices, nsource_colors, ntarget_colors)
    state = GraphComb.ColoredPortState(
        edges,
        zeros(Int, nvertices, nsource_colors),
        zeros(Int, nvertices, ntarget_colors),
    )
    canonical, _ = GraphComb.canonical_relabeling(problem, state)
    tuples = [
        (edge.source, edge.target, edge.source_color, edge.target_color) for
        edge in GraphComb.port_edges(canonical)
    ]
    sort!(tuples)
    return Tuple(tuples)
end

function brute_fermion_port_weights(
    problem,
    source_vertices,
    source_colors,
    target_vertices,
    target_colors,
    nvertices,
    nsource_colors,
    ntarget_colors,
)
    E = length(source_vertices)
    weights = Dict{Any,BigInt}()
    for permutation in permutations(1:E)
        edges = GraphComb.ColoredPortEdge[
            GraphComb.ColoredPortEdge(
                source_vertices[k],
                target_vertices[permutation[k]],
                source_colors[k],
                target_colors[permutation[k]],
            ) for k in 1:E
        ]
        key = canonical_port_edge_key(
            problem, edges, nvertices, nsource_colors, ntarget_colors
        )
        weight = BigInt(KC.pairing_sign(Fermion, permutation))
        weights[key] = get(weights, key, big(0)) + weight
    end
    filter!(pair -> !iszero(last(pair)), weights)
    return weights
end

function gc_fermion_port_weights(
    problem,
    source_vertices,
    source_colors,
    target_vertices,
    target_colors,
    nvertices,
    nsource_colors,
    ntarget_colors,
)
    transport = KC._GCFermionPortTransport(
        source_vertices, source_colors, target_vertices, target_colors
    )
    completions = GraphComb.generate_weighted(problem; transport)
    weights = Dict{Any,BigInt}()
    for completion in completions
        key = canonical_port_edge_key(
            problem, completion.edges, nvertices, nsource_colors, ntarget_colors
        )
        weights[key] = get(weights, key, big(0)) + completion.weight
    end
    filter!(pair -> !iszero(last(pair)), weights)
    return weights
end

function certify_fermion_port_problem(
    source_indices, target_indices, nvertices, nsource_colors, ntarget_colors
)
    source_vertices, source_colors = port_cells(source_indices, nsource_colors)
    target_vertices, target_colors = port_cells(target_indices, ntarget_colors)
    source_ports = port_counts(
        source_vertices, source_colors, nvertices, nsource_colors
    )
    target_ports = port_counts(
        target_vertices, target_colors, nvertices, ntarget_colors
    )
    compatibility = trues(nvertices, nsource_colors, nvertices, ntarget_colors)
    problem = GraphComb.ColoredPortProblem(
        ones(Int, nvertices), source_ports, target_ports, compatibility
    )

    expected = brute_fermion_port_weights(
        problem,
        source_vertices,
        source_colors,
        target_vertices,
        target_colors,
        nvertices,
        nsource_colors,
        ntarget_colors,
    )
    actual = gc_fermion_port_weights(
        problem,
        source_vertices,
        source_colors,
        target_vertices,
        target_colors,
        nvertices,
        nsource_colors,
        ntarget_colors,
    )
    return actual == expected
end

@testset "GraphCombinations weighted Wick backend" begin
    @testset "bosonic exactness" begin
        @qfields gcwick_ϕ::Boson
        c, q = gcwick_ϕ[Classical], gcwick_ϕ[Quantum]
        elastic = -(
            0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2)
        )

        expression = c(KC.Out()) * bar(q)(KC.In()) * elastic
        for term in KC.terms(expression)
            direct = KC._wick_contraction(
                term.args_nc, Val(3), Val(0); regularise=true, simplify=false
            )
            gc, stats = KC._gc_wick_contraction_with_stats(
                term.args_nc, Val(3), Val(0); regularise=true, simplify=false
            )
            @test signed_pairing_weights(gc) == signed_pairing_weights(direct)
            @test stats.canonicalization_calls >= 1
        end

        L = InteractionLagrangian(elastic)
        second_order = c(KC.Out()) * bar(c)(KC.In()) * L(1).lagrangian * L(2).lagrangian
        for term in KC.terms(second_order)
            direct = KC._wick_contraction(
                term.args_nc, Val(5), Val(1); regularise=true, simplify=true
            )
            gc, _ = KC._gc_wick_contraction_with_stats(
                term.args_nc, Val(5), Val(1); regularise=true, simplify=true
            )
            @test signed_pairing_weights(gc) == signed_pairing_weights(direct)
        end
    end

    @testset "fermionic cofactor cancellation" begin
        @qfields gcwick_pair_ψ::Fermion
        ψ₁ = gcwick_pair_ψ[One](KC.Bulk(1))
        ψ₂ = gcwick_pair_ψ[Two](KC.Bulk(1))
        bψ₂a = bar(gcwick_pair_ψ[Two])(KC.Bulk(2))
        bψ₂b = bar(gcwick_pair_ψ[Two])(KC.Bulk(3))
        term = ψ₁ * ψ₂ * bψ₂a * bψ₂b
        args = copy(KC.fields(term))

        direct = KC._wick_contraction(
            args, Val(2), Val(0); regularise=false, simplify=false
        )
        gc, stats = KC._gc_wick_contraction_with_stats(
            args, Val(2), Val(0); regularise=false, simplify=false
        )
        @test signed_pairing_weights(gc) == signed_pairing_weights(direct)
        @test isempty(gc)
        @test stats.automorphisms > 1
    end

    @testset "fermionic second-order relabeling transport" begin
        @qfields gcwick_ψ::Fermion
        ψ₁, ψ₂ = gcwick_ψ[One], gcwick_ψ[Two]
        vertex = ψ₁ * ψ₂ * bar(ψ₁) * bar(ψ₂)
        L = InteractionLagrangian(vertex)
        external_fields = KC.propagator_fields(L, nothing)
        products = KC.propagator_external_products(Fermion, external_fields...)

        saw_nontrivial_automorphism = false
        saw_merged_transition = false
        for in_out in products
            expression = in_out * L(1).lagrangian * L(2).lagrangian
            for term in KC.terms(expression)
                direct = KC._wick_contraction(
                    term.args_nc, Val(5), Val(1); regularise=false, simplify=false
                )
                gc, stats = KC._gc_wick_contraction_with_stats(
                    term.args_nc, Val(5), Val(1); regularise=false, simplify=false
                )
                @test signed_pairing_weights(gc) == signed_pairing_weights(direct)
                saw_nontrivial_automorphism |= stats.automorphisms > 1
                saw_merged_transition |= stats.merged_transitions > 0
            end
        end
        @test saw_nontrivial_automorphism
        @test saw_merged_transition
    end

    @testset "exhaustive fermionic port orientation" begin
        problem_count = 0
        for E in 1:3
            assignments = collect(Iterators.product(ntuple(_ -> 1:2, E)...))
            for source_indices in assignments, target_indices in assignments
                @test certify_fermion_port_problem(
                    source_indices, target_indices, 2, 1, 1
                )
                problem_count += 1
            end
        end
        for E in 1:2
            assignments = collect(Iterators.product(ntuple(_ -> 1:4, E)...))
            for source_indices in assignments, target_indices in assignments
                @test certify_fermion_port_problem(
                    source_indices, target_indices, 2, 2, 2
                )
                problem_count += 1
            end
        end
        @test problem_count == 356
    end
end
