# # Bosonic scattering and two-body loss
#
# Consider one complex bosonic field with coherent contact scattering and Markovian
# two-body loss. The point of this example is twofold: first obtain the physical
# collision kernel through the high-level API, then unpack exactly how the compiler
# derives it.

using KeldyshContraction
using KeldyshContraction: Regularisation

@qfields ϕ::Boson
c, q = ϕ[Classical], ϕ[Quantum]

# ## Interaction action
#
# For the coherent contact interaction we use the standard bosonic Keldysh vertex
# (up to the overall coupling `g`).

elastic = -(
    (1 // 2) * (c^2 + q^2) * bar(c) * bar(q) + (1 // 2) * c * q * (bar(c)^2 + bar(q)^2)
)

# The dissipative vertex must retain its finite Trotter ordering. The `+` and `-`
# labels distinguish the two sides of the equal-time contraction before the causal
# reduction removes that regulator structurally.

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
Lg = InteractionLagrangian(elastic, g)
Lγ = InteractionLagrangian(loss, γ)
L = Lg + Lγ

# `L` is the full microscopic interaction and keeps the `g`, `γ`, `g^2`, `gγ`, ...
# sectors distinct. When only one sector is required, it is cheaper to build the
# corresponding propagator from the relevant process directly; no unrelated diagrams
# then need to be generated.
#
# ## Ask for the collision kernel
#
# A user who wants the kinetic result does not need to drive the compiler manually.
# Build the desired perturbative propagator and call `collision_kernel`.

Gγ = DressedPropagator(Lγ, Val(1), Val(3); preserve_regularisation=true)
Cγ = collision_kernel(Gγ)
Cγ

# In the normalization used here the generated first-order loss kernel is
#
# ```math
# C_n^{(\gamma)}(k)=-4\gamma\int_q n_k n_q.
# ```
#
# The compact call above is the ordinary API. The rest of this section exposes the
# same computation stage by stage. Each returned value has a physics-facing LaTeX
# display, so the rendered documentation shows the represented equation rather than
# an internal field dump.
#
# ## 1. Fourier transform
#
# Coordinate derivatives become exact momentum polynomials and the two-point diagrams
# acquire routed frequency-momentum variables.

GF = fourier_transform(Gγ)
GF

# ## 2. One-particle-irreducible self-energy
#
# `SelfEnergy` extracts the fixed-order 1PI contribution and amputates the external
# propagators. No self-consistent Dyson iteration is implied.

ΣF = SelfEnergy(GF)
ΣF

# ## 3. Wigner representation
#
# Centre/relative coordinates are converted to Wigner variables. The present kinetic
# compiler uses the homogeneous zeroth-gradient collision term.

ΣW = wigner_transform(ΣF; gradient_order=Val(0))
ΣW

# ## 4. Spectral/statistical kinetic representation
#
# The Keldysh propagators are expressed through spectral information and the
# statistical distribution, `G^K=-iFA` at zeroth gradient order.

KΣ = kinetic_expression(ΣW)
KΣ

# ## 5. Complete Kadanoff--Baym collision identity
#
# The off-shell collision expression is assembled before any on-shell reduction,
#
# ```math
# I_{\rm coll}=i\Sigma^K-F_k A_\Sigma,
# \qquad
# A_\Sigma=i(\Sigma^R-\Sigma^A).
# ```

I = off_shell_collision_expression(KΣ)
I

# ## 6. Spectral/dispersive decomposition
#
# Retarded and advanced lines are decomposed as
# `G^R=D-iA/2` and `G^A=D+iA/2`. This keeps both shell and principal-value physics.

SD = spectral_dispersive_collision(I)
SD

# ## 7. Causal frequency reduction
#
# Internal frequencies are eliminated exactly. Finite sectors retain
# `δ(ΔE)` and, when present, `PV(1/ΔE)`. Genuine pinches remain explicit blockers.

R = reduce_frequency_collision(SD)
R

# ## 8. Occupation reduction
#
# For bosons, `F=1+2n`; after the complete Kadanoff--Baym combination has been formed,
# the statistical polynomial is converted to the Bose gain/loss polynomial.

N = occupation_reduced_expression(R)
N

# ## 9. Loop-momentum quotient
#
# Dummy loop variables related by signed permutations represent the same physical
# integral. They are quotiented only after the occupation polynomial is known.

Q = quotient_loop_momenta(N)
Q

# ## 10. Final physical kernel
#
# The explicit compiler route terminates in the same `CollisionKernel` as the one-line
# user API above.

Cγ_explicit = collision_kernel(Q)
Cγ_explicit

@assert KeldyshContraction.collision_kernel_terms(Cγ) ==
    KeldyshContraction.collision_kernel_terms(Cγ_explicit)

# ## Elastic two-body scattering
#
# The identical compiler applies at second order in the coherent coupling. Two quartic
# vertices plus the two external legs give five propagator edges at this order.

Gg2 = DressedPropagator(Lg, Val(2), Val(5); preserve_regularisation=true)
Cg2 = collision_kernel(Gg2)
Cg2

# Its Bose-enhanced gain/loss structure is
#
# ```math
# C_n^{(g^2)}(k)\propto
# 2\!\left[(1+n_k)(1+n_p)n_qn_r
# -n_kn_p(1+n_q)(1+n_r)\right]\delta(\Delta E),
# ```
#
# with the exact momentum-routing and phase-space factors retained by the symbolic
# kernel. In the full multi-process expansion, mixed `g\gamma` sectors can contain
# principal-value contributions. If a selected sector also contains a genuine pinch,
# `collision_kernel(G)` refuses to hide it: follow the explicit route through
# `R = reduce_frequency_collision(...)` and inspect `reduced_blocked_terms(R)` instead.
