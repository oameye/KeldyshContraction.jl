include("gc_wick_adapter.jl")

using Statistics

# Causal code for an allowed port pair:
# 0 = noncausal, 1 = source vertex -> target vertex (advanced),
# 2 = target vertex -> source vertex (retarded).
struct FastCausalPortPruner
    causal_codes::Array{UInt8,4}
end

function FastCausalPortPruner(problem, lookup)
    codes = zeros(UInt8, size(problem.compatibility))
    for (cell, contraction) in lookup
        propagator = KC.propagator_type(contraction...)
        if KC.is_advanced(propagator)
            codes[cell...] = 0x01
        elseif KC.is_retarded(propagator)
            codes[cell...] = 0x02
        end
    end
    return FastCausalPortPruner(codes)
end

@inline function fast_causal_code(pruner::FastCausalPortPruner, edge::GC._PortEdge)::UInt8
    return pruner.causal_codes[
        edge.source, edge.source_color, edge.target, edge.target_color
    ]
end

@inline function fast_causal_endpoints(edge::GC._PortEdge, code::UInt8)::Tuple{Int,Int}
    return code == 0x01 ? (edge.source, edge.target) : (edge.target, edge.source)
end

function has_fast_causal_path(
    pruner::FastCausalPortPruner,
    edges::Vector{GC._PortEdge},
    current::Int,
    target::Int,
    last_index::Int,
    visited::UInt128,
)::Bool
    current == target && return true
    current > 128 && error("KC causal pruning currently assumes at most 128 graph vertices")
    current_bit = one(UInt128) << (current - 1)
    !iszero(visited & current_bit) && return false
    visited |= current_bit

    @inbounds for i in 1:last_index
        edge = edges[i]
        code = fast_causal_code(pruner, edge)
        iszero(code) && continue
        source, destination = fast_causal_endpoints(edge, code)
        source == destination && continue
        source == current || continue
        destination == target && return true
        has_fast_causal_path(pruner, edges, destination, target, last_index, visited) &&
            return true
    end
    return false
end

function (pruner::FastCausalPortPruner)(state::GC._PortMatchingState)::Bool
    isempty(state.edges) && return true
    new_index = length(state.edges)
    new_edge = state.edges[new_index]
    code = fast_causal_code(pruner, new_edge)
    iszero(code) && return true

    source, destination = fast_causal_endpoints(new_edge, code)
    if source == destination
        @inbounds for i in 1:(new_index - 1)
            edge = state.edges[i]
            old_code = fast_causal_code(pruner, edge)
            iszero(old_code) && continue
            old_source, old_destination = fast_causal_endpoints(edge, old_code)
            old_source == source && old_destination == source && return false
        end
        return true
    end

    return !has_fast_causal_path(
        pruner, state.edges, destination, source, new_index - 1, zero(UInt128)
    )
end

struct FastCausalPruningStats
    terms::Int
    transitions::Int
    canonicalizations::Int
    pruned_transitions::Int
    completed_states::Int
    filter_survivors::Int
    canonical_unique::Int
end
FastCausalPruningStats() = FastCausalPruningStats(0, 0, 0, 0, 0, 0, 0)
Base.:+(a::FastCausalPruningStats, b::FastCausalPruningStats) = FastCausalPruningStats(
    a.terms + b.terms,
    a.transitions + b.transitions,
    a.canonicalizations + b.canonicalizations,
    a.pruned_transitions + b.pruned_transitions,
    a.completed_states + b.completed_states,
    a.filter_survivors + b.filter_survivors,
    a.canonical_unique + b.canonical_unique,
)

function fast_causal_pruned_adapter_weights(
    args_nc::Vector{KC.Field{S}},
    ::Val{E};
    regularise=true,
    _set_reg_to_zero=false,
    simplify=false,
) where {S<:KC.Statistics,E}
    problem, lookup = build_gc_problem(args_nc, Val(E); regularise, _set_reg_to_zero)
    pruner = FastCausalPortPruner(problem, lookup)
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

    stats = FastCausalPruningStats(
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

function fast_pruning_terms(
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

function validate_fast_causal_component(
    in_out,
    L,
    ::Val{O},
    ::Val{E};
    simplify=true,
    _set_reg_to_zero=true,
) where {O,E}
    terms = fast_pruning_terms(in_out, L, Val(O))
    regularise = KC.should_regularise(L.lagrangian)
    total = FastCausalPruningStats()
    for args_nc in terms
        expected, _ = adapter_weights(
            args_nc,
            Val(E);
            regularise,
            _set_reg_to_zero,
            simplify,
        )
        actual, stats = fast_causal_pruned_adapter_weights(
            args_nc,
            Val(E);
            regularise,
            _set_reg_to_zero,
            simplify,
        )
        actual == expected || error("incremental causal pruning changed exact KC weights")
        total += stats
    end
    return total
end

function print_fast_stats(workload, component, stats)
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

function validate_fast_causal(workload, L, order, edges)
    products = KC.propagator_external_products(Boson, KC.propagator_fields(L, nothing)...)
    total = FastCausalPruningStats()
    for (name, in_out) in zip(("K", "R", "A"), products)
        stats = validate_fast_causal_component(in_out, L, Val(order), Val(edges))
        print_fast_stats(workload, name, stats)
        total += stats
    end
    print_fast_stats(workload, "total", total)
    return total
end

function timed_fast_samples(f, samples::Int)
    result = Vector{NamedTuple{(:time, :bytes),Tuple{Float64,Int}}}(undef, samples)
    for i in 1:samples
        measurement = @timed f()
        result[i] = (time=measurement.time, bytes=measurement.bytes)
    end
    return result
end

function print_fast_timing(workload, component, backend, outputs, samples)
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

function benchmark_fast_causal_component(
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
    terms = fast_pruning_terms(in_out, L, Val(O))
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
            weights, _ = fast_causal_pruned_adapter_weights(
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
    unpruned_outputs == pruned_outputs || error("incremental causal pruning changed output count")

    print_fast_timing(
        workload,
        component,
        "gc_unpruned",
        unpruned_outputs,
        timed_fast_samples(unpruned, samples),
    )
    print_fast_timing(
        workload,
        component,
        "gc_fast_causal_pruned",
        pruned_outputs,
        timed_fast_samples(pruned, samples),
    )
    return nothing
end

function benchmark_fast_causal(workload, L, order, edges; samples=3)
    products = KC.propagator_external_products(Boson, KC.propagator_fields(L, nothing)...)
    for (name, in_out) in zip(("K", "R", "A"), products)
        benchmark_fast_causal_component(
            workload, name, in_out, L, Val(order), Val(edges); samples
        )
    end
    return nothing
end

println("# incremental causal-pruning acceptance and search statistics")
println(
    "workload\tcomponent\tterms\ttransitions\tcanonicalizations\tpruned_transitions\t",
    "completed_states\tfilter_survivors\tcanonical_unique",
)
validate_fast_causal("boson_g3", L_g, 3, 7)
validate_fast_causal("boson_g4", L_g, 4, 9)

println("# incremental causal-pruned GC timing")
println("workload\tcomponent\tbackend\toutputs\tbest_s\tmedian_s\tbest_bytes\tmedian_bytes")
benchmark_fast_causal("boson_g3", L_g, 3, 7; samples=3)
benchmark_fast_causal("boson_g4", L_g, 4, 9; samples=3)
