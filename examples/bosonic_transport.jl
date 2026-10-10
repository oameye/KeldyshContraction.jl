# # Bosonic scattering and two-body loss
#
# Consider one complex bosonic field with coherent contact scattering and Markovian
# two-body loss. We first ask for the physical collision kernel directly, then follow
# the same calculation through the objects that carry the derivation.

using KeldyshContraction
using KeldyshContraction: Regularisation

@qfields ϕ::Boson
c, q = ϕ[Classical], ϕ[Quantum]

# ## Microscopic interaction
#
# For the coherent contact interaction we use the standard bosonic Keldysh vertex,
# with the coupling kept outside the field polynomial.

elastic = -(
    (1 // 2) * (c^2 + q^2) * bar(c) * bar(q) + (1 // 2) * c * q * (bar(c)^2 + bar(q)^2)
)

# The dissipative vertex retains its finite Trotter ordering. The `+` and `-` labels
# distinguish the two sides of an equal-time contraction until the causal reduction
# removes the regulator structurally.

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

# The full action `L` keeps the perturbative sectors `g`, `γ`, `g^2`, `gγ`, ...
# distinct. If only one sector is needed, constructing from the corresponding process
# avoids generating unrelated diagrams.
#
# ## The physical loss kernel
#
# For the first-order loss problem the ordinary user-facing calculation is just the
# perturbative propagator followed by `collision_kernel`.

Gγ = DressedPropagator(Lγ, Val(1), Val(3); preserve_regularisation=true)
Cγ = collision_kernel(Gγ)
Cγ

# In the normalization used here the result is
#
# ```math
# C_n^{(\gamma)}(k)=-4\gamma\int_q n_k n_q.
# ```
#
# The compact call is the endpoint. To see where that kernel comes from, we now expose
# the same derivation without changing the calculation.
#
# ## From diagrams to the kinetic equation
#
# Fourier transformation assigns exact momentum variables to the coordinate-space
# contractions. The resulting object already shows the routed propagator expression.

GF = fourier_transform(Gγ)
GF

# The 1PI self-energy is then obtained by amputating the external propagators and
# removing reducible two-point diagrams. This is a fixed-order perturbative
# self-energy, not a self-consistent Dyson solution.

ΣF = SelfEnergy(GF)
ΣF

# Passing to Wigner variables separates centre and relative coordinates. For the
# homogeneous collision problem we keep the zeroth-gradient contribution.

ΣW = wigner_transform(ΣF; gradient_order=Val(0))
ΣW

# The kinetic representation replaces internal Keldysh propagators by spectral and
# statistical objects. At this order,
#
# ```math
# G^K=-iFA.
# ```

KΣ = kinetic_expression(ΣW)
KΣ

# The full Kadanoff--Baym collision combination is formed before any shell projection,
#
# ```math
# I_{\rm coll}=i\Sigma^K-F_k A_\Sigma,
# \qquad
# A_\Sigma=i(\Sigma^R-\Sigma^A).
# ```
#
# This ordering matters: cancellations between the self-energy components belong to
# the collision identity itself and must occur before frequency reduction.

I = off_shell_collision_expression(KΣ)
I

# ## Resolving the causal structure
#
# Retarded and advanced propagators are decomposed into spectral and dispersive parts,
#
# ```math
# G^R=D-\frac{i}{2}A,
# \qquad
# G^A=D+\frac{i}{2}A.
# ```
#
# The displayed object below is the resulting concrete `A/D/F` integrand, not merely a
# label for this transformation.

SD = spectral_dispersive_collision(I)
SD

# Internal frequencies can now be integrated exactly. Regular terms retain their
# physical shell or principal-value support, while genuine pinches remain explicit
# blockers rather than being assigned a finite value.

R = reduce_frequency_collision(SD)
R

# Only after this causal reduction do we replace the bosonic statistical function by
# occupations, `F=1+2n`. The displayed polynomial is therefore already the physical
# Bose gain/loss structure carried by each reduced sector.

N = occupation_reduced_expression(R)
N

# Dummy loop variables related by signed permutations describe the same momentum
# integral. Quotienting them here gives a canonical representative without changing
# the physical occupation or kinematic factors.

Q = quotient_loop_momenta(N)
Q

# Lowering the canonical sectors produces the same physical kernel as the direct call
# with which we started.

Cγ_explicit = collision_kernel(Q)
Cγ_explicit

@assert KeldyshContraction.collision_kernel_terms(Cγ) ==
    KeldyshContraction.collision_kernel_terms(Cγ_explicit)

# ## Coherent two-body scattering
#
# The same compiler applies at second order in the coherent coupling. Two quartic
# vertices plus the external two-point legs give five propagator edges.

Gg2 = DressedPropagator(Lg, Val(2), Val(5); preserve_regularisation=true)
Cg2 = collision_kernel(Gg2)
Cg2

# The result has the expected Bose-enhanced gain/loss structure,
#
# ```math
# C_n^{(g^2)}(k)\propto
# 2\!\left[(1+n_k)(1+n_p)n_qn_r
# -n_kn_p(1+n_q)(1+n_r)\right]\delta(\Delta E),
# ```
#
# with the exact routing and phase-space factors retained by the symbolic kernel. In a
# full multi-process expansion, mixed `g\gamma` sectors may also carry principal-value
# contributions. If a selected sector contains a genuine pinch,
# `collision_kernel(G)` refuses to hide it: stop at `R`, inspect
# `reduced_regular_terms(R)` and `reduced_blocked_terms(R)`, and make that approximation
# choice explicitly.
