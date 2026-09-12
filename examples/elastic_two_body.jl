# ## Elastic Two Body Scattering
using KeldyshContraction

# ## System

# The interaction action of elastic two-body scattering is
# ```math
# S_\mathrm{int} = -\frac{g}{2} \int d^d x \,
# [ (\bar{\phi}_+\phi_+)^2 - (\bar{\phi}_-\phi_-)^2 ].
# ```
# This is the standard local interaction for s-wave bosonic scattering. In the RAK basis,
# ```math
# S_\mathrm{int} = -\frac{g}{2} \int d^d x \,
# [\bar{\phi}_c\bar{\phi}_q(\phi_c^2+\phi_q^2)
# +\phi_c\phi_q(\bar{\phi}_c^2+\bar{\phi}_q^2)].
# ```

@qfields ϕ::Boson
c, q = ϕ[Classical], ϕ[Quantum]

elastic2boson = -(
    (1 // 2) * (c^2 + q^2) * bar(c) * bar(q) +
    (1 // 2) * c * q * (bar(c)^2 + bar(q)^2)
)
L_int = InteractionLagrangian(elastic2boson, :g)

# The normalization identity ``Z=1`` is enforced by the Wick-contraction machinery and is
# covered by the package's diagrammatic tests. The supported user workflow starts from the
# interaction and constructs the requested propagator correction with static perturbation
# order and edge count.

# ## First order

G1 = DressedPropagator(L_int, Val(1), Val(3))
Σ1 = SelfEnergy(G1)

# For this elastic interaction the first-order self-energy has no collision contribution:
# ``Σ^K=0`` and ``Σ^R=Σ^A``. The nontrivial elastic Boltzmann kernel appears at second order.

# ## Second-order kinetic kernel

G2 = DressedPropagator(L_int, Val(2), Val(5))
GF2 = fourier_transform(G2)
ΣF2 = SelfEnergy(GF2)
ΣW2 = wigner_transform(ΣF2; gradient_order=Val(0))
kinetic2 = kinetic_expression(ΣW2)
off_shell2 = off_shell_collision_expression(kinetic2)
spectral2 = spectral_dispersive_collision(off_shell2)
reduced2 = reduce_frequency_collision(spectral2)
occupation2 = occupation_reduced_expression(reduced2)
quotiented2 = quotient_loop_momenta(occupation2)
kernel2 = collision_kernel(quotiented2)

# `kernel2` is the canonical on-shell collision kernel after exact frequency reduction,
# occupation conversion, and physical loop-momentum quotienting.
kernel2
