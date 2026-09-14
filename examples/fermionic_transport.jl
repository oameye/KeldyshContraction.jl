# # Fermionic p-wave scattering and two-body loss
#
# For identical spinless fermions the local s-wave channel vanishes.  A minimal
# odd-wave process therefore carries one spatial derivative.  We use one Cartesian
# channel,
#
# ```math
# P_x=\psi\,\partial_x\psi,
# ```
#
# and combine coherent p-wave scattering with two-body loss.

using KeldyshContraction
using KeldyshContraction: Regularisation

@qfields ψ::Fermion
ψ1, ψ2 = ψ[One], ψ[Two]

# ## Interaction action
#
# Fermions use the asymmetric Larkin--Ovchinnikov rotation.  The barred fields are
# independent Grassmann variables, and their branch reconstruction differs from the
# unbarred one.

sqrt2 = sqrt(2.0)
ψplus = (ψ1 + ψ2) / sqrt2
ψminus = (ψ1 - ψ2) / sqrt2
barψplus = (bar(ψ2) + bar(ψ1)) / sqrt2
barψminus = (bar(ψ2) - bar(ψ1)) / sqrt2

Pplus = ψplus * partial(ψplus, :x)
Pminus = ψminus * partial(ψminus, :x)
barPplus = partial(barψplus, :x) * barψplus
barPminus = partial(barψminus, :x) * barψminus

# The coherent vertex is the contour difference of `P_x† P_x`.

elastic = -(1 // 2) * (barPplus * Pplus - barPminus * Pminus)

# The loss vertex keeps the finite-Trotter ordering inherited from the microscopic
# Lindblad step.  The regulator is not an optional numerical convention: it fixes the
# equal-time causal contractions before the strict kinetic reduction.

plus = Regularisation.Plus
minus = Regularisation.Minus
Pplus_prev = ψplus(minus) * partial(ψplus(minus), :x)
Pminus_prev = ψminus(minus) * partial(ψminus(minus), :x)
barPplus_next = partial(barψplus(plus), :x) * barψplus(plus)
barPminus_next = partial(barψminus(plus), :x) * barψminus(plus)
loss =
    (1 // 2) *
    im *
    (
        (barPplus_next - barPminus_next) * (Pplus_prev + Pminus_prev) +
        (barPplus_next + barPminus_next) * (Pplus_prev - Pminus_prev)
    )

@syms gp γp
L = InteractionLagrangian(elastic, gp) + InteractionLagrangian(loss, γp)

# ## Ask for the collision kernel
#
# As for bosons, the ordinary API is one call once the perturbative sector is selected.

G1 = DressedPropagator(L, Val(1), Val(3); preserve_regularisation=true)
Gγp = G1[γp]
Cγp = collision_kernel(Gγp)
Cγp

# The generated first-order odd-wave loss kernel is
#
# ```math
# C_n^{(\gamma_p)}(k)
# =-\gamma_p\int_q (k_x-q_x)^2 n_k n_q.
# ```
#
# The relative-momentum factor is produced by the derivative vertices and survives all
# subsequent canonical reductions.
#
# ## Inspecting the derivation
#
# We now unpack the same `collision_kernel(Gγp)` call.  The stages are identical to the
# bosonic compiler; statistics enters through Grassmann signs and the final map
# `F_F=1-2n`.
#
# ### 1. Fourier transform

GF = fourier_transform(Gγp)
GF

# Derivatives become exact odd-wave momentum polynomials while Wick permutation parity
# has already fixed the fermionic signs.
#
# ### 2. One-particle-irreducible self-energy

ΣF = SelfEnergy(GF)
ΣF

# The result is a fixed-order 1PI self-energy.  No nonlinear Dyson iteration is hidden
# in this constructor.
#
# ### 3. Wigner representation

ΣW = wigner_transform(ΣF; gradient_order=Val(0))
ΣW

# We retain the homogeneous zeroth-gradient collision side of the kinetic equation.
#
# ### 4. Spectral/statistical kinetic representation

KΣ = kinetic_expression(ΣW)
KΣ

# At this stage internal Keldysh lines are parameterized by spectral functions and the
# fermionic statistical distribution.
#
# ### 5. Complete Kadanoff--Baym collision identity

I = off_shell_collision_expression(KΣ)
I

# The compiler forms
#
# ```math
# I_{\rm coll}=i\Sigma^K-F_kA_\Sigma,
# \qquad A_\Sigma=i(\Sigma^R-\Sigma^A),
# ```
#
# before any shell projection.  This ordering is essential for the cancellation and
# survival of the correct causal sectors.
#
# ### 6. Spectral/dispersive decomposition

SD = spectral_dispersive_collision(I)
SD

# Retarded and advanced propagators are decomposed into spectral and dispersive pieces;
# principal-value information is preserved.
#
# ### 7. Causal frequency reduction

R = reduce_frequency_collision(SD)
R

# The first-order loss sector is regular.  At higher orders the same object also reports
# genuine causal pinches explicitly instead of assigning them a finite value.
#
# ### 8. Occupation reduction

N = occupation_reduced_expression(R)
N

# For fermions `F=1-2n`; Pauli blocking therefore emerges here from the same generic
# statistical compiler used for bosons.
#
# ### 9. Loop-momentum quotient

Q = quotient_loop_momenta(N)
Q

# Dummy loop labels are canonically quotiented only after the occupation polynomial and
# exact odd-wave kinematic factor are known.
#
# ### 10. Final physical kernel

Cγp_explicit = collision_kernel(Q)
Cγp_explicit

@assert KeldyshContraction.collision_kernel_terms(Cγp) ==
    KeldyshContraction.collision_kernel_terms(Cγp_explicit)

# ## Second-order p-wave scattering
#
# The regular coherent `g_p^2` sector can again be compiled directly.

G2 = DressedPropagator(L, Val(2), Val(3); preserve_regularisation=true)
Cgp2 = collision_kernel(G2[gp ^ 2])
Cgp2

# The certified shell contribution is
#
# ```math
# C_n^{(g_p^2)}
# =\frac{g_p^2}{8}(k_x-p_x)^2(q_x-r_x)^2
# \Big[(1-n_k)(1-n_p)n_qn_r
# -n_kn_p(1-n_q)(1-n_r)\Big]\delta(\Delta E).
# ```
#
# Canonical frequency and loop reduction produces nine physical kernel sectors and no
# causal blocker in this channel.
#
# ## Dissipative second order and mixed sectors
#
# The pure-loss shell contribution is
#
# ```math
# C_n^{(\gamma_p^2)}
# =\frac{\gamma_p^2}{8}(k_x-p_x)^2(q_x-r_x)^2
# \Big[(1-n_k)(1-n_p)n_qn_r
# +n_kn_p(1-n_q)(1-n_r)
# -2\,n_k n_p n_q n_r\Big]\delta(\Delta E).
# ```
#
# It has nine finite shell sectors together with 84 genuine causal-frequency pinches.
# Consequently `collision_kernel(G2[γp^2])` is intentionally not a silent "finite-part"
# operation: the high-level compiler refuses to discard those blockers.  To study this
# sector, follow the explicit pipeline through `R = reduce_frequency_collision(...)`,
# inspect `reduced_regular_terms(R)` and `reduced_blocked_terms(R)`, and only then choose
# whether a finite regular branch is the intended approximation.
#
# The mixed coherent--dissipative contribution is purely dispersive in the certified
# strict reduction,
#
# ```math
# C_n^{(g_p\gamma_p)}
# =-\frac{g_p\gamma_p}{2}(k_x-p_x)^2(q_x-r_x)^2
# n_kn_p(1-n_q-n_r)\operatorname{PV}\!\left(\frac{1}{\Delta E}\right).
# ```
#
# There is no mixed shell/Born term: the finite branch contains nine PV sectors, while
# 30 genuine pinches remain explicit blockers.
