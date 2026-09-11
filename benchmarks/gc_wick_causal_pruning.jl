include("gc_wick_adapter.jl")

using Statistics

mutable struct CausalPortPruner{S<:KC.Statistics}
    lookup::Dict{NTuple{4,Int},KC.Contraction{S}}
    contractions::Vector{KC.Contraction{S}}
end

function (pruner::CausalPortPruner{S})(state::GC._PortMatchingState)::Bool where {S}
    resize!(pruner.contractions, length(state.edges))
    for i in eachindex(state.edges)
        edge = state.edges[i]
        pruner.contractions[i] = pruner.lookup[
            (edge.source, edge.source_color, edge.target, edge.target_color)
        ]
    end
    return !KC.has_zero_loop(pruner.contractions)
end

struct CausalPruningStats
    terms::Int
    transitions::Int
    canonicalizations::Int
    pruned_transitions::Int
    completed_states::Int
    filter_survivors::Int
    canonical_unique::Int
end
CausalPruningStats() = CausalPruningStats(0, 0, 0, 0, 0, 0, 0)
Base.:+(a::CausalPruningStats, b::CausalPruningStats) = CausalPruningStats(
    a.terms + b.terms,
    a.transitions + b.transitions,
    a.canonicalizations + b.canonicalizations,
    a.pruned_transitions + b.pruned_transitions,
    a.completed_states + b.completed_states,
    a.filter_survivors + b.filter_survivors,
    a.canonical_unique + b.canonical_unique,
)

function causal_pruned_adapter_weights(
    args_nc::Vector{KC.Field{S}},
    ::Val{E};
    regularise=true,
    _set_reg_to_zero=false,
    simplify=false,
) where {S<:KC.Statistics,E}
    problem, lookup = build_gc_problem(args_nc, Val(E); regularise, _set_reg_to_zero)
    pruner = CausalPortPruner(lookup, KC.Contraction{S}[])
    results, search = GC._weighted_port_matchings_pruned_with_stats(problem, pruner)

    weights = Dict{FixedVector{E,KC.Contraction{S}},BigInt}()
    filter_survivors = 0
    for (edges, multiplicity) in results
        contractions = Vector{KC.Contraction{S}}(undef, E)
        for i in eachindex(edges)
            edge = edges[i]
            contractions[i] = lookup[
                (edge.source, edge.source_color, edge.target, edge.target_color)
            ]
        end
        KC.passes_wick_filters(contractions) || continue
        filter_survivors += 1
        final_contractions, sign = simplify ?
            KC.advanced_to_retarded(contractions, 1) : (contractions, 1)
        key = KC.sorted_wick_key(KC.canonicalize(final_contractions), Val(E))
        weights[key] = get(weights, key, big(0)) + multiplicity * Int(sign)
    end
    filter!(pair -> !iszero(last(pair)), weights)

    stats = CausalPruningStats(
        1,
        search.transitions,
        search.canonicalization_calls,
        search.pruned_transitions,
        length(results),
        filter_survivors,
        length(weights),
    )
    return weights, stats
end

function pruning_workload_terms(
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

function validate_causal_pruning_component(
    in_out,
    L,
    ::Val{O},
    ::Val{E};
    simplify=true,
    _set_reg_to_zero=true,
) where {O,E}
    terms = pruning_workload_terms(in_out, L, Val(O))
    regularise = KC.should_regularise(L.lagrangian)
    total = CausalPruningStats()
    for args_nc in terms
        expected, _ = adapter_weights(
            args_nc,
            Val(E);
            regularise,
            _set_reg_to_zero,
            simplify,
        )
        actual, stats = causal_pruned_adapter_weights(
            args_nc,
            Val(E);
            regularise,
            _set_reg_to_zero,
            simplify,
        )
        actual == expected || error("causal-pruned GC adapter changed exact KC weights")
        total += stats
    end
    return total
end

function print_pruning_stats(workload, component, stats)
    println(
        join(
            (
                workload,
                component,
                stats.terms,
                stats.transitions,
                stats.canonicalizations,
                stats.pruned_transitions,
                stats.completed_states,
                stats.filter_survivors,
                stats.canonical_unique,
            ),
            '\t',
        ),
    )
end

function validate_causal_pruning(workload, L, order, edges)
    products = KC.propagator_external_products(Boson, KC.propagator_fields(L, nothing)...)
    total = CausalPruningStats()
    for (name, in_out) in zip(("K", "R", "A"), products)
        stats = validate_causal_pruning_component(in_out, L, Val(order), Val(edges))
        print_pruning_stats(workload, name, stats)
        total += stats
    end
    print_pruning_stats(workload, "total", total)
    return total
end

function timed_pruning_samples(f, samples::Int)
    result = Vector{NamedTuple{(:time, :bytes),Tuple{Float64,Int}}}(undef, samples)
    for i in 1:samples
        measurement = @timed f()
        result[i] = (time=measurement.time, bytes=measurement.bytes)
    end
    return result
end

function print_pruning_timing(workload, component, backend, outputs, samples)
    times = getproperty.(samples, :time)
    bytes = getproperty.(samples, :bytes)
    println(
        join(
            (
                workload,
                component,
                backend,
                outputs,
                minimum(times),
                median(times),
                minimum(bytes),
                median(bytes),
            ),
            '\t',
        ),
    )
end

function benchmark_causal_pruning_component(
    workload,
    component,
    in_out,
    L,
    ::Val{O},
    ::Val{E};
    samples=3,
    simplify=true,
    _set_reg_to_zero=true,
) where {O,E}
    terms = pruning_workload_terms(in_out, L, Val(O))
    regularise = KC.should_regularise(L.lagrangian)

    unpruned = () -> begin
        outputs = 0
        for args_nc in terms
            weights, _ = adapter_weights(
                args_nc,
                Val(E);
                regularise,
                _set_reg_to_zero,
                simplify,
            )
            outputs += length(weights)
        end
        outputs
    end
    pruned = () -> begin
        outputs = 0
        for args_nc in terms
            weights, _ = causal_pruned_adapter_weights(
                args_nc,
                Val(E);
                regularise,
                _set_reg_to_zero,
                simplify,
            )
            outputs += length(weights)
        end
        outputs
    end

    unpruned_outputs = unpruned()
    pruned_outputs = pruned()
    unpruned_outputs == pruned_outputs || error("causal pruning changed output count")

    print_pruning_timing(
        workload,
        component,
        "gc_unpruned",
        unpruned_outputs,
        timed_pruning_samples(unpruned, samples),
    )
    print_pruning_timing(
        workload,
        component,
        "gc_causal_pruned",
        pruned_outputs,
        timed_pruning_samples(pruned, samples),
    )
    return nothing
end

function benchmark_causal_pruning(workload, L, order, edges; samples=3)
    products = KC.propagator_external_products(Boson, KC.propagator_fields(L, nothing)...)
    for (name, in_out) in zip(("K", "R", "A"), products)
        benchmark_causal_pruning_component(
            workload, name, in_out, L, Val(order), Val(edges); samples
        )
    end
    return nothing
end

println("# exact causal-pruning acceptance and search statistics")
println(
    "workload\tcomponent\tterms\ttransitions\tcanonicalizations\tpruned_transitions\t",
    "completed_states\tfilter_survivors\tcanonical_unique",
)
validate_causal_pruning("boson_g3", L_g, 3, 7)
validate_causal_pruning("boson_g4", L_g, 4, 9)

println("# causal-pruned GC timing")
println("workload\tcomponent\tbackend\toutputs\tbest_s\tmedian_s\tbest_bytes\tmedian_bytes")
benchmark_causal_pruning("boson_g3", L_g, 3, 7; samples=3)
benchmark_causal_pruning("boson_g4", L_g, 4, 9; samples=3)
