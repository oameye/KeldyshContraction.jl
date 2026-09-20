using KeldyshContraction
import KeldyshContraction as KC

function self_energy_exact(a, b)
    return a.keldysh == b.keldysh &&
           a.retarded == b.retarded &&
           a.advanced == b.advanced &&
           KC.parameters(a) == KC.parameters(b) &&
           KC.target_family(a) == KC.target_family(b) &&
           KC.order(a) == KC.order(b)
end

function best_self_energy_measurement(f; samples=2)
    f()
    best_time = Inf
    best_bytes = typemax(Int)
    for _ in 1:samples
        measurement = @timed f()
        best_time = min(best_time, measurement.time)
        best_bytes = min(best_bytes, measurement.bytes)
    end
    return best_time, best_bytes
end

function benchmark_self_energy(name, L, ::Val{O}, ::Val{E}; samples=2) where {O,E}
    baseline = () -> SelfEnergy(DressedPropagator(L, Val(O), Val(E)))
    direct = () -> KC._direct_self_energy(L, Val(O), Val(E))

    reference = baseline()
    candidate = direct()
    self_energy_exact(candidate, reference) ||
        error("direct self-energy semantic mismatch: $name")

    baseline_time, baseline_bytes = best_self_energy_measurement(baseline; samples)
    direct_time, direct_bytes = best_self_energy_measurement(direct; samples)

    outputs =
        length(candidate.keldysh) + length(candidate.retarded) + length(candidate.advanced)
    println(
        join(
            (
                name,
                O,
                E,
                outputs,
                baseline_time,
                direct_time,
                direct_time / baseline_time,
                baseline_bytes,
                direct_bytes,
                direct_bytes / baseline_bytes,
            ),
            '\t',
        ),
    )
    return nothing
end

@qfields direct_se_ϕ::Boson
c, q = direct_se_ϕ[Classical], direct_se_ϕ[Quantum]
elastic = -(0.5 * (c^2 + q^2) * bar(c) * bar(q) + 0.5 * c * q * (bar(c)^2 + bar(q)^2))
L_g = InteractionLagrangian(elastic, :g)

loss =
    0.5 *
    bar(c) *
    bar(q) *
    (
        c(KC.Regularisation.Minus) * c(KC.Regularisation.Minus) +
        q(KC.Regularisation.Minus) * q(KC.Regularisation.Minus)
    ) -
    0.5 *
    c(KC.Regularisation.Plus) *
    q(KC.Regularisation.Plus) *
    (bar(c) * bar(c) + bar(q) * bar(q)) +
    bar(c) *
    bar(q) *
    (
        c(KC.Regularisation.Plus) * q(KC.Regularisation.Plus) +
        c(KC.Regularisation.Minus) * q(KC.Regularisation.Minus)
    )
L_γ = InteractionLagrangian(loss, :γ)

println(
    "workload\torder\tedges\toutputs\tbaseline_s\tdirect_s\ttime_ratio\tbaseline_bytes\tdirect_bytes\talloc_ratio",
)
benchmark_self_energy("boson_g2", L_g, Val(2), Val(5); samples=3)
benchmark_self_energy("boson_g3", L_g, Val(3), Val(7); samples=2)
benchmark_self_energy("boson_g4", L_g, Val(4), Val(9); samples=1)
benchmark_self_energy("boson_gamma2", L_γ, Val(2), Val(5); samples=3)
benchmark_self_energy("boson_gamma3", L_γ, Val(3), Val(7); samples=2)
