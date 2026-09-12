# ## Two-Body Loss
using KeldyshContraction
using KeldyshContraction: Regularisation

# ## System and equal-time regularisation

# For local two-body loss with jump ``L=\sqrt{\gamma}\,\phi^2``, the dissipative interaction
# action requires the finite Trotter shifts inherited from the discrete-time coherent-state
# construction. In the RAK basis we keep those shifts explicitly until the frequency reduction
# has used their equal-time information.

@qfields ϕ::Boson
c, q = ϕ[Classical], ϕ[Quantum]
plus = Regularisation.Plus
minus = Regularisation.Minus

loss2boson =
    (1 // 2) *
    im *
    (
        bar(c) * bar(q) * (c(minus)^2 + q(minus)^2) -
        c(plus) * q(plus) * (bar(c)^2 + bar(q)^2) +
        2 * bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
    )

L_int = InteractionLagrangian(loss2boson, :γ)

# `preserve_regularisation=true` is the public control for equal-time/Trotter-sensitive
# calculations. The default high-level propagator construction removes these shifts after the
# Wick contraction, which is appropriate when no later reduction needs their provenance.

# ## End-to-end first-order loss kernel

G = DressedPropagator(
    L_int, Val(1), Val(3); simplify=true, preserve_regularisation=true
)
GF = fourier_transform(G)
ΣF = SelfEnergy(GF)
ΣW = wigner_transform(ΣF; gradient_order=Val(0))
kinetic = kinetic_expression(ΣW)
off_shell = off_shell_collision_expression(kinetic)
spectral = spectral_dispersive_collision(off_shell)
reduced = reduce_frequency_collision(spectral)
occupation = occupation_reduced_expression(reduced)
quotiented = quotient_loop_momenta(occupation)
kernel = collision_kernel(quotiented)

# The resulting occupation kernel is the generated two-body-loss law. In the package's
# normalization its first-order integrand is ``-4γ n_k n_q``.
kernel
