include("gc_wick_adapter.jl")

function split_workload_terms(
    in_out::KC.QMul{CI,Boson},
    L::KC.InteractionLagrangian{CL,Boson},
    ::Val{O},
) where {CI<:Number,CL<:Number,O}
    result = Vector{Vector{KC.Field{Boson}}}()
    for coefficients in Combinatorics.multiexponents(length(L.lagrangian), O)
        idxs = KC.indices_from_counts(coefficients)
        mult = Combinatorics.multinomial(coefficients...)
        qmul = mult * prod(L(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))
        push!(result, (in_out * qmul).args_nc)
    end
    return result
end

function materialize_contractions(
    edges,
    lookup::Dict{NTuple{4,Int},KC.Contraction{S}},
    ::Val{E},
) where {S<:KC.Statistics,E}
    contractions = Vector{KC.Contraction{S}}(undef, E)
    for i in eachindex(edges)
        edge = edges[i]
        contractions[i] = lookup[
            (edge.source, edge.source_color, edge.target, edge.target_color)
        ]
    end
    return contractions
end

function filter_split(
    terms,
    ::Val{E};
    regularise=true,
    _set_reg_to_zero=true,
) where {E}
    total = 0
    disconnected_only = 0
    causal_only = 0
    both = 0
    survivors = 0

    for args_nc in terms
        problem, lookup = build_gc_problem(
            args_nc, Val(E); regularise, _set_reg_to_zero
        )
        results, _ = GC._weighted_port_matchings_with_stats(problem)
        for (edges, _) in results
            contractions = materialize_contractions(edges, lookup, Val(E))
            connected = KC.is_connected(contractions)
            causal_zero = KC.has_zero_loop(contractions)
            total += 1
            if connected && !causal_zero
                survivors += 1
            elseif !connected && causal_zero
                both += 1
            elseif !connected
                disconnected_only += 1
            else
                causal_only += 1
            end
        end
    end

    total == disconnected_only + causal_only + both + survivors ||
        error("filter split does not partition completed quotient states")
    return (; total, disconnected_only, causal_only, both, survivors)
end

function print_filter_split(workload, component, split)
    println(
        join(
            (
                workload,
                component,
                split.total,
                split.disconnected_only,
                split.causal_only,
                split.both,
                split.survivors,
            ),
            '\t',
        ),
    )
end

function measure_filter_split(workload, L, order, edges)
    products = KC.propagator_external_products(Boson, KC.propagator_fields(L, nothing)...)
    for (name, in_out) in zip(("K", "R", "A"), products)
        terms = split_workload_terms(in_out, L, Val(order))
        regularise = KC.should_regularise(L.lagrangian)
        split = filter_split(
            terms, Val(edges); regularise, _set_reg_to_zero=true
        )
        print_filter_split(workload, name, split)
    end
    return nothing
end

println("# completed quotient-state filter split")
println("workload\tcomponent\ttotal\tdisconnected_only\tcausal_only\tboth\tsurvivors")
measure_filter_split("boson_g3", L_g, 3, 7)
measure_filter_split("boson_g4", L_g, 4, 9)
