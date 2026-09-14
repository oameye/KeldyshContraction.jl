"""
    collision_kernel(result::ReducedFrequencyCollision)

Lower a completely regular reduced collision to the final `CollisionKernel`.

This method is deliberately strict. Genuine causal-frequency blockers and unresolved
finite-Trotter states do not have a finite value in the strict quasiparticle reduction and
are therefore never discarded implicitly. Inspect them with `reduced_blocked_terms` and
`reduced_trotter_terms`. If only the finite regular branch is wanted explicitly, call
`collision_kernel(occupation_reduced_expression(result))`.
"""
function collision_kernel(result::ReducedFrequencyCollision)
    blocked = reduced_blocked_terms(result)
    unresolved = reduced_trotter_terms(result)
    if !isempty(blocked) || !isempty(unresolved)
        throw(
            ArgumentError(
                "cannot construct a complete CollisionKernel: the reduced collision contains " *
                "$(length(blocked)) causal-frequency blocker(s) and $(length(unresolved)) " *
                "unresolved finite-Trotter state(s). Inspect reduced_blocked_terms(result) and " *
                "reduced_trotter_terms(result), or explicitly lower only the finite regular " *
                "branch with collision_kernel(occupation_reduced_expression(result)).",
            ),
        )
    end
    return collision_kernel(occupation_reduced_expression(result))
end

"""
    collision_kernel(G::DressedPropagator; gradient_order=Val(0))

Compile one perturbative dressed-propagator sector directly to its physical strict-QP collision
kernel.

This is the ordinary user entry point for the kinetic compiler. Internally it performs the
certified sequence

`Fourier -> 1PI self-energy -> Wigner -> kinetic R/A/K -> Kadanoff--Baym collision ->`
` spectral/dispersive decomposition -> causal frequency reduction -> F-to-n reduction ->`
` dummy-loop quotient -> CollisionKernel`.

The individual transformations remain public so the derivation can be inspected stage by stage.
The high-level route is intentionally equivalent to that explicit pipeline. It also inherits the
strict blocker policy of `collision_kernel(::ReducedFrequencyCollision)`: genuine pinches or
unresolved finite-Trotter states are reported rather than silently dropped.

Only gradient order zero is currently implemented by the Wigner/kinetic stack.
"""
function collision_kernel(G::DressedPropagator; gradient_order::Val=Val(0))
    GF = fourier_transform(G)
    ΣF = SelfEnergy(GF)
    ΣW = wigner_transform(ΣF; gradient_order=gradient_order)
    KΣ = kinetic_expression(ΣW)
    I = off_shell_collision_expression(KΣ)
    SD = spectral_dispersive_collision(I)
    R = reduce_frequency_collision(SD)
    return collision_kernel(R)
end
