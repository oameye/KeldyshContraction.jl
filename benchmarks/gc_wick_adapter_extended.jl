include("gc_wick_adapter.jl")

function validate_mixed_component(
    in_out::KC.QMul{CI,Boson},
    L1::KC.InteractionLagrangian{C1,Boson},
    L2::KC.InteractionLagrangian{C2,Boson};
    simplify=true,
    _set_reg_to_zero=true,
) where {CI<:Number,C1<:Number,C2<:Number}
    total = AdapterStats()
    qadd = 2 * L1(1).lagrangian * L2(2).lagrangian
    regularise = KC.should_regularise(qadd)
    for arg in KC.terms(qadd)
        total += validate_fields(
            (in_out * arg).args_nc,
            Val(5),
            Val(KC.max_edges(2));
            regularise,
            _set_reg_to_zero,
            simplify,
        )
    end
    return total
end

function validate_mixed(workload, L1, L2; simplify=true, _set_reg_to_zero=true)
    products = KC.propagator_external_products(Boson, KC.propagator_fields(L1, nothing)...)
    total = AdapterStats()
    for (name, in_out) in zip(("K", "R", "A"), products)
        stats = validate_mixed_component(
            in_out, L1, L2; simplify, _set_reg_to_zero
        )
        print_stats(workload, name, stats)
        total += stats
    end
    print_stats(workload, "total", total)
    return total
end

println("# extended bosonic acceptance")
validate_interaction("boson_gamma1", L_γ, 1, 3; simplify=true, _set_reg_to_zero=true)
validate_mixed("boson_g_gamma", L_g, L_γ; simplify=true, _set_reg_to_zero=true)
validate_interaction("boson_g4", L_g, 4, 9)
