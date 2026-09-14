# # Bosonic scattering and two-body loss
#
# This example starts from a single Schwinger--Keldysh action containing coherent contact
# scattering and Markovian two-body loss and follows both processes to kinetic collision
# kernels.

using KeldyshContraction
using KeldyshContraction: Regularisation

@qfields ϕ::Boson
c, q = ϕ[Classical], ϕ[Quantum]

# ## Interaction action
#
# For contact scattering, `H_int = (g/2) ϕ†²ϕ²` gives the RAK interaction

elastic =
    -(1 // 2) * ((c^2 + q^2) * bar(c) * bar(q) + c * q * (bar(c)^2 + bar(q)^2))

# Two-body loss requires the finite-Trotter ordering of the jump vertex. The regularisation
# labels are kept until the structural equal-time reduction.

plus = Regularisation.Plus
minus = Regularisation.Minus
loss =
    (1 // 2) *
    im *
    (
        bar(c) * bar(q) * (c(minus)^2 + q(minus)^2) -
        c(plus) * q(plus) * (bar(c)^2 + bar(q)^2) +
        2 * bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
    )

@syms g γ
L = InteractionLagrangian(elastic, g) + InteractionLagrangian(loss, γ)

# ## Self-energy to collision kernel
#
# The same public pipeline is used for every perturbative sector.

function kinetic_kernel(G)
    GF = fourier_transform(G)
    ΣF = SelfEnergy(GF)
    ΣW = wigner_transform(ΣF; gradient_order=Val(0))
    KΣ = kinetic_expression(ΣW)
    I = off_shell_collision_expression(KΣ)
    SD = spectral_dispersive_collision(I)
    R = reduce_frequency_collision(SD)
    N = occupation_reduced_expression(R)
    Q = quotient_loop_momenta(N)
    return collision_kernel(Q)
end

# The first-order dissipative sector gives the loss kernel.

G1 = DressedPropagator(L, Val(1), Val(3); preserve_regularisation=true)
Cγ = kinetic_kernel(G1[γ])

# With the convention used in this action,
# ```math
# C_n^{(\gamma)}(k)=-4\gamma\int_q n_k n_q.
# ```
Cγ

# The coherent collision integral starts at second order.

G2 = DressedPropagator(L, Val(2), Val(5); preserve_regularisation=true)
Cg2 = kinetic_kernel(G2[g^2])

# Its on-shell occupation factor is the standard bosonic gain--loss combination,
# ```math
# 2\big[(1+n_k)(1+n_p)n_qn_r
#      -n_kn_p(1+n_q)(1+n_r)\big]\,\delta(\Delta E),
# ```
# with momentum conservation fixing the fourth momentum. The same second-order expansion also
# retains the mixed `gγ` principal-value sector rather than forcing it onto shell.
Cg2
