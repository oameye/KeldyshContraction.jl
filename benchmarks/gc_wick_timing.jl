include("gc_wick_adapter.jl")

using Statistics

function workload_terms(
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

function old_matcher(
    terms,
    ::Val{E},
    ::Val{E2};
    regularise=true,
    _set_reg_to_zero=true,
    simplify=true,
) where {E,E2}
    outputs = 0
    for args_nc in terms
        outputs += length(
            KC._wick_contraction(
                args_nc,
                Val(E),
                Val(E2);
                regularise,
                _set_reg_to_zero,
                simplify,
            ),
        )
    end
    return outputs
end

function gc_matcher(
    terms,
    ::Val{E};
    regularise=true,
    _set_reg_to_zero=true,
    simplify=true,
) where {E}
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
    return outputs
end

function gc_build_stage(
    terms,
    ::Val{E};
    regularise=true,
    _set_reg_to_zero=true,
) where {E}
    return [
        build_gc_problem(
            args_nc, Val(E); regularise, _set_reg_to_zero
        ) for args_nc in terms
    ]
end

function gc_generate_stage(prepared)
    return [GC._weighted_port_matchings_with_stats(first(input))[1] for input in prepared]
end

function gc_postprocess_results(
    results,
    lookup::Dict{NTuple{4,Int},KC.Contraction{S}},
    ::Val{E};
    simplify=true,
) where {S<:KC.Statistics,E}
    weights = Dict{FixedVector{E,KC.Contraction{S}},BigInt}()
    for (edges, multiplicity) in results
        contractions = Vector{KC.Contraction{S}}(undef, E)
        for i in eachindex(edges)
            edge = edges[i]
            contractions[i] = lookup[
                (edge.source, edge.source_color, edge.target, edge.target_color)
            ]
        end
        KC.passes_wick_filters(contractions) || continue
        final_contractions, sign = simplify ?
            KC.advanced_to_retarded(contractions, 1) : (contractions, 1)
        key = KC.sorted_wick_key(KC.canonicalize(final_contractions), Val(E))
        weights[key] = get(weights, key, big(0)) + multiplicity * Int(sign)
    end
    filter!(pair -> !iszero(last(pair)), weights)
    return length(weights)
end

function gc_postprocess_stage(prepared, generated, ::Val{E}; simplify=true) where {E}
    outputs = 0
    for i in eachindex(prepared, generated)
        outputs += gc_postprocess_results(
            generated[i], last(prepared[i]), Val(E); simplify
        )
    end
    return outputs
end

function timed_samples(f, samples::Int)
    measurements = Vector{NamedTuple{(:time, :bytes),Tuple{Float64,Int}}}(undef, samples)
    for i in 1:samples
        result = @timed f()
        measurements[i] = (time=result.time, bytes=result.bytes)
    end
    return measurements
end

function print_timing(workload, component, backend, outputs, samples)
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

function benchmark_component(
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
    terms = workload_terms(in_out, L, Val(O))
    regularise = KC.should_regularise(L.lagrangian)
    old = () -> old_matcher(
        terms,
        Val(E),
        Val(KC.max_edges(O));
        regularise,
        _set_reg_to_zero,
        simplify,
    )
    new = () -> gc_matcher(
        terms,
        Val(E);
        regularise,
        _set_reg_to_zero,
        simplify,
    )

    # Compile both paths before collecting measurements. The exact semantic oracle for these
    # workloads is exercised separately by gc_wick_adapter_extended.jl on the same pinned SHAs.
    old_outputs = old()
    new_outputs = new()
    old_outputs == new_outputs || error("timed matcher output counts disagree")

    old_samples = timed_samples(old, samples)
    new_samples = timed_samples(new, samples)
    print_timing(workload, component, "old", old_outputs, old_samples)
    print_timing(workload, component, "gc", new_outputs, new_samples)
    return nothing
end

function benchmark_gc_stages(
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
    terms = workload_terms(in_out, L, Val(O))
    regularise = KC.should_regularise(L.lagrangian)
    prepared = gc_build_stage(
        terms, Val(E); regularise, _set_reg_to_zero
    )
    generated = gc_generate_stage(prepared)
    expected_outputs = gc_matcher(
        terms,
        Val(E);
        regularise,
        _set_reg_to_zero,
        simplify,
    )
    gc_postprocess_stage(prepared, generated, Val(E); simplify) == expected_outputs ||
        error("GC stage decomposition changed matcher output counts")

    build = () -> gc_build_stage(
        terms, Val(E); regularise, _set_reg_to_zero
    )
    generate = () -> gc_generate_stage(prepared)
    postprocess = () -> gc_postprocess_stage(prepared, generated, Val(E); simplify)

    build()
    generate()
    postprocess()
    print_timing(workload, component, "gc_build", length(terms), timed_samples(build, samples))
    print_timing(
        workload,
        component,
        "gc_generate",
        length(generated),
        timed_samples(generate, samples),
    )
    print_timing(
        workload,
        component,
        "gc_postprocess",
        expected_outputs,
        timed_samples(postprocess, samples),
    )
    return nothing
end

function benchmark_interaction(workload, L, order, edges; samples=3)
    products = KC.propagator_external_products(Boson, KC.propagator_fields(L, nothing)...)
    for (name, in_out) in zip(("K", "R", "A"), products)
        benchmark_component(
            workload, name, in_out, L, Val(order), Val(edges); samples
        )
    end
    return nothing
end

function benchmark_interaction_stages(workload, L, order, edges; samples=3)
    products = KC.propagator_external_products(Boson, KC.propagator_fields(L, nothing)...)
    for (name, in_out) in zip(("K", "R", "A"), products)
        benchmark_gc_stages(
            workload, name, in_out, L, Val(order), Val(edges); samples
        )
    end
    return nothing
end

println("# matcher-only warmed timings")
println("workload\tcomponent\tbackend\toutputs\tbest_s\tmedian_s\tbest_bytes\tmedian_bytes")
benchmark_interaction("boson_g3", L_g, 3, 7; samples=3)
benchmark_interaction("boson_g4", L_g, 4, 9; samples=3)

println("# GC stage decomposition")
println("workload\tcomponent\tbackend\toutputs\tbest_s\tmedian_s\tbest_bytes\tmedian_bytes")
benchmark_interaction_stages("boson_g3", L_g, 3, 7; samples=3)
benchmark_interaction_stages("boson_g4", L_g, 4, 9; samples=3)
