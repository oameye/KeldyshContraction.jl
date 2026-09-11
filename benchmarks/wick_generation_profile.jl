using KeldyshContraction
using Combinatorics
using SmallCollections: FixedVector

import KeldyshContraction as KC

mutable struct MatchingStats
    terms::Int
    candidate_rows::Int
    candidate_slots::Int
    max_candidate_width::Int
    recursive_states::Int
    complete_matchings::Int
    raw_unique::Int
    raw_zero_weights::Int
    filter_survivors::Int
    canonicalization_calls::Int
    canonical_unique::Int
    canonical_zero_weights::Int
end

MatchingStats() = MatchingStats(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

function add!(a::MatchingStats, b::MatchingStats)
    a.terms += b.terms
    a.candidate_rows += b.candidate_rows
    a.candidate_slots += b.candidate_slots
    a.max_candidate_width = max(a.max_candidate_width, b.max_candidate_width)
    a.recursive_states += b.recursive_states
    a.complete_matchings += b.complete_matchings
    a.raw_unique += b.raw_unique
    a.raw_zero_weights += b.raw_zero_weights
    a.filter_survivors += b.filter_survivors
    a.canonicalization_calls += b.canonicalization_calls
    a.canonical_unique += b.canonical_unique
    a.canonical_zero_weights += b.canonical_zero_weights
    return a
end

function profile_matching!(
    f,
    candidates::Tuple,
    contractions::Vector,
    permutation::Vector{Int},
    used::Vector{Bool},
    k::Int,
    stats::MatchingStats,
)
    stats.recursive_states += 1
    if k > length(candidates)
        stats.complete_matchings += 1
        f(contractions, permutation)
        return nothing
    end

    for (l, contraction) in candidates[k]
        used[l] && continue
        used[l] = true
        contractions[k] = contraction
        permutation[k] = l
        profile_matching!(f, candidates, contractions, permutation, used, k + 1, stats)
        used[l] = false
    end
    return nothing
end

function profile_fields(
    args_nc::Vector{KC.Field{S}},
    ::Val{E},
    ::Val{E2};
    regularise=true,
    _set_reg_to_zero=false,
    simplify=false,
) where {S<:KC.Statistics,E,E2}
    stats = MatchingStats()
    stats.terms = 1

    destroys, creates = KC.prepare_args(args_nc, Val(E))
    ps = map(KC.position, args_nc)
    skip = KC.has_in(ps) && KC.has_out(ps)
    candidates = KC.wick_candidates(
        destroys,
        creates,
        Val(E);
        regularise,
        _set_reg_to_zero,
        skip_external_pair=skip,
    )

    widths = map(length, candidates)
    stats.candidate_rows = E
    stats.candidate_slots = sum(widths)
    stats.max_candidate_width = isempty(widths) ? 0 : maximum(widths)

    candidate_vector_type = fieldtype(typeof(candidates), 1)
    candidate_type = fieldtype(eltype(candidate_vector_type), 2)
    contractions = Vector{candidate_type}(undef, E)
    permutation = Vector{Int}(undef, E)
    used = fill(false, E)

    matching_weights = Dict{NTuple{E,UInt8},Int}()
    profile_matching!(candidates, contractions, permutation, used, 1, stats) do cs, perm
        key = KC.wick_matching_key(cs, candidates, Val(E))
        weight = Int(KC.pairing_sign(S, perm))
        matching_weights[key] = get(matching_weights, key, 0) + weight
        return nothing
    end

    stats.raw_unique = length(matching_weights)
    stats.raw_zero_weights = count(iszero, values(matching_weights))

    canonical_weights = Dict{FixedVector{E,KC.Contraction{S}},Int}()
    for (key, weight) in matching_weights
        iszero(weight) && continue
        cs = KC.contractions_from_matching_key(key, candidates, KC.Contraction{S})
        KC.passes_wick_filters(cs) || continue
        stats.filter_survivors += 1

        final_cs, simplification_sign = if simplify
            KC.advanced_to_retarded(cs, 1)
        else
            cs, 1
        end
        stats.canonicalization_calls += 1
        canonical = KC.canonicalize(final_cs)
        canonical_key = KC.sorted_wick_key(canonical, Val(E))
        final_weight = weight * Int(simplification_sign)
        canonical_weights[canonical_key] =
            get(canonical_weights, canonical_key, 0) + final_weight
    end

    stats.canonical_unique = count(x -> !iszero(x), values(canonical_weights))
    stats.canonical_zero_weights = count(iszero, values(canonical_weights))
    return stats
end

function profile_component(
    in_out::KC.QMul{CI,S},
    L::KC.InteractionLagrangian{CL,S},
    ::Val{O},
    ::Val{E};
    simplify=true,
    _set_reg_to_zero=true,
) where {CI<:Number,CL<:Number,S<:KC.Statistics,O,E}
    total = MatchingStats()
    regularise = KC.should_regularise(L.lagrangian)
    l = length(L.lagrangian)

    for coefficients in Combinatorics.multiexponents(l, O)
        idxs = KC.indices_from_counts(coefficients)
        mult = Combinatorics.multinomial(coefficients...)
        qmul = mult * prod(L(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))
        term = in_out * qmul
        add!(
            total,
            profile_fields(
                term.args_nc,
                Val(E),
                Val(KC.max_edges(O));
                regularise,
                _set_reg_to_zero,
                simplify,
            ),
        )
    end
    return total
end

function profile_mixed_component(
    in_out::KC.QMul{CI,S},
    L1::KC.InteractionLagrangian{C1,S},
    L2::KC.InteractionLagrangian{C2,S};
    simplify=true,
    _set_reg_to_zero=true,
) where {CI<:Number,C1<:Number,C2<:Number,S<:KC.Statistics}
    total = MatchingStats()
    qadd = 2 * L1(1).lagrangian * L2(2).lagrangian
    regularise = KC.should_regularise(qadd)
    for arg in KC.terms(qadd)
        term = in_out * arg
        add!(
            total,
            profile_fields(
                term.args_nc,
                Val(5),
                Val(KC.max_edges(2));
                regularise,
                _set_reg_to_zero,
                simplify,
            ),
        )
    end
    return total
end

function statistics_type(::KC.InteractionLagrangian{C,S}) where {C,S}
    return S
end

function print_header()
    println(
        "workload\tcomponent\tterms\tcandidate_rows\tcandidate_slots\tmax_width\t",
        "recursive_states\tcomplete_matchings\traw_unique\traw_zero\tfilter_survivors\t",
        "canonicalization_calls\tcanonical_unique\tcanonical_zero",
    )
end

function print_row(workload, component, s::MatchingStats)
    println(
        join(
            (
                workload,
                component,
                s.terms,
                s.candidate_rows,
                s.candidate_slots,
                s.max_candidate_width,
                s.recursive_states,
                s.complete_matchings,
                s.raw_unique,
                s.raw_zero_weights,
                s.filter_survivors,
                s.canonicalization_calls,
                s.canonical_unique,
                s.canonical_zero_weights,
            ),
            "\t",
        ),
    )
end

function profile_interaction(
    workload,
    L::KC.InteractionLagrangian,
    order::Int,
    edges::Int;
    simplify=true,
    _set_reg_to_zero=true,
)
    S = statistics_type(L)
    fields = KC.propagator_fields(L, nothing)
    products = KC.propagator_external_products(S, fields...)
    total = MatchingStats()
    for (name, in_out) in zip(("K", "R", "A"), products)
        stats = profile_component(
            in_out,
            L,
            Val(order),
            Val(edges);
            simplify,
            _set_reg_to_zero,
        )
        print_row(workload, name, stats)
        add!(total, stats)
    end
    print_row(workload, "total", total)
    return total
end

function profile_mixed(workload, L1, L2; simplify=true, _set_reg_to_zero=true)
    S = statistics_type(L1)
    fields = KC.propagator_fields(L1, nothing)
    products = KC.propagator_external_products(S, fields...)
    total = MatchingStats()
    for (name, in_out) in zip(("K", "R", "A"), products)
        stats = profile_mixed_component(
            in_out, L1, L2; simplify, _set_reg_to_zero
        )
        print_row(workload, name, stats)
        add!(total, stats)
    end
    print_row(workload, "total", total)
    return total
end

function timed_dressed(label, L, order, edges; kwargs...)
    f = () -> DressedPropagator(L, Val(order), Val(edges); kwargs...)
    f() # compile/warm up
    samples = [@timed f() for _ in 1:3]
    times = [sample.time for sample in samples]
    best = samples[argmin(times)]
    println(
        "TIMING\t",
        label,
        "\tseconds=",
        best.time,
        "\tbytes=",
        best.bytes,
        "\tfinal_diagrams=",
        sum(length, (best.value.keldysh, best.value.retarded, best.value.advanced)),
    )
    return nothing
end

print_header()

@qfields ϕprofile::Boson
c, q = ϕprofile[Classical], ϕprofile[Quantum]

elastic = -(
    0.5 * (c^2 + q^2) * bar(c) * bar(q) +
    0.5 * c * q * (bar(c)^2 + bar(q)^2)
)
L_g = InteractionLagrangian(elastic, :g)

loss =
    0.5 * bar(c) * bar(q) *
    (c(KC.Regularisation.Minus) * c(KC.Regularisation.Minus) +
     q(KC.Regularisation.Minus) * q(KC.Regularisation.Minus)) -
    0.5 * c(KC.Regularisation.Plus) * q(KC.Regularisation.Plus) *
    (bar(c) * bar(c) + bar(q) * bar(q)) +
    bar(c) * bar(q) *
    (c(KC.Regularisation.Plus) * q(KC.Regularisation.Plus) +
     c(KC.Regularisation.Minus) * q(KC.Regularisation.Minus))
L_γ = InteractionLagrangian(loss, :γ)

profile_interaction("boson_g2", L_g, 2, 5)
profile_interaction("boson_gamma1", L_γ, 1, 3; simplify=true, _set_reg_to_zero=true)
profile_interaction("boson_gamma2", L_γ, 2, 5; simplify=true, _set_reg_to_zero=true)
profile_mixed("boson_g_gamma", L_g, L_γ; simplify=true, _set_reg_to_zero=true)
profile_interaction("boson_g3", L_g, 3, 7)

timed_dressed("boson_g2", L_g, 2, 5)
timed_dressed("boson_gamma2", L_γ, 2, 5; simplify=true, _set_reg_to_zero=true)
timed_dressed("boson_g3", L_g, 3, 7)

@qfields ψprofile::Fermion
ψ₁, ψ₂ = ψprofile[One], ψprofile[Two]
fermion_vertex = ψ₁ * ψ₂ * bar(ψ₁) * bar(ψ₂)
L_f = InteractionLagrangian(fermion_vertex, :u)
profile_interaction("fermion_quartic2", L_f, 2, 5; simplify=false)
timed_dressed("fermion_quartic2", L_f, 2, 5; simplify=false)

∂xψ₂ = partial(ψ₂, :x)
derivative_vertex = ψ₁ * ∂xψ₂ * bar(ψ₁) * bar(∂xψ₂)
L_p = InteractionLagrangian(derivative_vertex, :γ)
profile_interaction("fermion_derivative1", L_p, 1, 3; simplify=false)
timed_dressed("fermion_derivative1", L_p, 1, 3; simplify=false)
