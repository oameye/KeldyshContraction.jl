# # Fermionic p-wave scattering and two-body loss
#
# For identical spinless fermions the local s-wave channel vanishes. A minimal
# odd-wave process therefore carries one spatial derivative. We use one Cartesian
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
bψ1, bψ2 = bar(ψ1), bar(ψ2)

# ## Microscopic p-wave processes
#
# Fermions use the asymmetric Larkin--Ovchinnikov inverse rotation. We keep the common
# `1/√2` factors outside the branch sums. Each pair contains two such factors, giving
# the exact coherent coefficient `-1/8`. Derivatives act on the fundamental fields
# before the branch sums are formed, making the linearity of `∂x` explicit in the
# symbolic representation.

ψplus = ψ1 + ψ2
∂ψplus = partial(ψ1, :x) + partial(ψ2, :x)
ψminus = ψ1 - ψ2
∂ψminus = partial(ψ1, :x) - partial(ψ2, :x)

bψplus = bψ1 + bψ2
∂bψplus = partial(bψ1, :x) + partial(bψ2, :x)
bψminus = bψ2 - bψ1
∂bψminus = partial(bψ2, :x) - partial(bψ1, :x)

Pplus = ψplus * ∂ψplus
Pminus = ψminus * ∂ψminus
Pplus_dagger = ∂bψplus * bψplus
Pminus_dagger = ∂bψminus * bψminus

elastic = -(1 // 8) * (Pplus_dagger * Pplus - Pminus_dagger * Pminus)

# The loss vertex uses the finite-Trotter contour ordering required by the dissipative
# theory. The incoming `+` branch is evaluated one step back, the forward branch
# difference carries the `+` regulator, and the barred pair is unshifted.

plus = Regularisation.Plus
minus = Regularisation.Minus

ψplus_minus = ψ1(minus) + ψ2(minus)
∂ψplus_minus = partial(ψ1(minus), :x) + partial(ψ2(minus), :x)
ψplus_plus = ψ1(plus) + ψ2(plus)
∂ψplus_plus = partial(ψ1(plus), :x) + partial(ψ2(plus), :x)
ψminus_plus = ψ1(plus) - ψ2(plus)
∂ψminus_plus = partial(ψ1(plus), :x) - partial(ψ2(plus), :x)

Pplus_minus = ψplus_minus * ∂ψplus_minus
Pplus_plus = ψplus_plus * ∂ψplus_plus
Pminus_plus = ψminus_plus * ∂ψminus_plus

loss =
    (1 // 8) *
    im *
    (
        (Pplus_dagger - Pminus_dagger) * Pplus_minus -
        (Pplus_plus - Pminus_plus) * Pminus_dagger
    )

@syms gp γp
Lg = InteractionLagrangian(elastic, gp)
Lγ = InteractionLagrangian(loss, γp)
L = Lg + Lγ

# The full interaction `L` keeps coherent, dissipative, and mixed perturbative sectors
# distinct. For a focused calculation we can compile one process directly and avoid
# generating unrelated sectors.
#
# ## The first-order loss kernel
#
# The physical first-order result is obtained from the perturbative propagator with one
# high-level compiler call.

Gγp = DressedPropagator(Lγ, Val(1), Val(3); preserve_regularisation=true)
Cγp = collision_kernel(Gγp)
Cγp

# The generated odd-wave loss kernel is
#
# ```math
# C_n^{(\gamma_p)}(k)
# =-\gamma_p\int_q (k_x-q_x)^2 n_k n_q.
# ```
#
# The relative-momentum factor is not inserted by hand. It comes from the derivative
# vertices and survives the diagrammatic, causal, statistical, and loop reductions.
# To see that explicitly, we now follow the same computation stage by stage.
#
# ## From derivative vertices to the kinetic equation
#
# Fourier transformation converts the spatial derivatives into exact odd-wave momentum
# polynomials. Wick permutation parity has already fixed the fermionic signs, so the
# rendered propagator contains both the routed lines and the derivative kinematics.

GF = fourier_transform(Gγp)
GF

# The fixed-order 1PI self-energy is then extracted by amputating the external lines and
# rejecting reducible two-point diagrams. There is no nonlinear Dyson iteration hidden
# in this constructor.

ΣF = SelfEnergy(GF)
ΣF

# Wigner transformation introduces centre and relative variables. We keep the
# homogeneous zeroth-gradient collision term used by the present kinetic compiler.

ΣW = wigner_transform(ΣF; gradient_order=Val(0))
ΣW

# Internal Keldysh lines are next expressed through the spectral function and the
# fermionic statistical distribution. At this order the relation is
#
# ```math
# G^K=-iFA.
# ```

KΣ = kinetic_expression(ΣW)
KΣ

# The collision side is assembled using the complete Kadanoff--Baym identity,
#
# ```math
# I_{\rm coll}=i\Sigma^K-F_kA_\Sigma,
# \qquad A_\Sigma=i(\Sigma^R-\Sigma^A).
# ```
#
# This combination is formed before any shell projection so that the correct
# cancellations and surviving causal sectors are determined by the field theory rather
# than by a later approximation.

I = off_shell_collision_expression(KΣ)
I

# ## From causal propagators to occupations
#
# Retarded and advanced propagators are separated into spectral and dispersive pieces,
#
# ```math
# G^R=D-\frac{i}{2}A,
# \qquad
# G^A=D+\frac{i}{2}A.
# ```
#
# The object below therefore displays the actual routed `A/D/F` collision expression,
# including the odd-wave momentum factors and finite-Trotter shifts.

SD = spectral_dispersive_collision(I)
SD

# Exact internal-frequency integration then produces the regular shell or
# principal-value sectors. The first-order loss channel is regular; at higher orders
# genuine causal pinches are retained explicitly instead of being assigned a finite
# value.

R = reduce_frequency_collision(SD)
R

# Only after the causal reduction do we replace the fermionic statistical function by
# occupations, `F=1-2n`. This is where the statistical polynomial becomes the physical
# Pauli gain/loss polynomial.

N = occupation_reduced_expression(R)
N

# Signed permutations of dummy loop variables represent the same integral. The quotient
# chooses a canonical loop basis while preserving the occupation polynomial and the
# exact p-wave kinematics.

Q = quotient_loop_momenta(N)
Q

# The canonical sectors finally lower to the same collision kernel returned by the
# high-level call at the start of the example.

Cγp_explicit = collision_kernel(Q)
Cγp_explicit

@assert KeldyshContraction.collision_kernel_terms(Cγp) ==
    KeldyshContraction.collision_kernel_terms(Cγp_explicit)

# ## Coherent p-wave scattering
#
# The regular coherent `g_p^2` sector is compiled by the same route. Two quartic
# vertices and the external two-point legs give five propagator edges.

Ggp2 = DressedPropagator(Lg, Val(2), Val(5); preserve_regularisation=true)
Cgp2 = collision_kernel(Ggp2)
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
# Expanding the full `L = Lg + Lγ` interaction to second order also produces
# `g_p\gamma_p` and `\gamma_p^2` sectors. Both contain genuine blocked causal support,
# so they are not appropriate for a silent one-line finite-part extraction.
#
# The pure-loss regular shell branch is
#
# ```math
# C_n^{(\gamma_p^2)}
# =\frac{\gamma_p^2}{8}(k_x-p_x)^2(q_x-r_x)^2
# \Big[(1-n_k)(1-n_p)n_qn_r
# +n_kn_p(1-n_q)(1-n_r)
# -2\,n_k n_p n_q n_r\Big]\delta(\Delta E).
# ```
#
# It contains nine finite shell sectors together with 84 genuine causal-frequency
# pinches. `collision_kernel` therefore refuses to discard the blocked branch. Stop at
# `R`, inspect `reduced_regular_terms(R)` and `reduced_blocked_terms(R)`, and decide
# explicitly whether retaining only the finite branch is the intended approximation.
#
# The mixed coherent--dissipative regular contribution is purely dispersive,
#
# ```math
# C_n^{(g_p\gamma_p)}
# =-\frac{g_p\gamma_p}{2}(k_x-p_x)^2(q_x-r_x)^2
# n_kn_p(1-n_q-n_r)\operatorname{PV}\!\left(\frac{1}{\Delta E}\right).
# ```
#
# There is no mixed shell/Born term: the finite branch contains nine principal-value
# sectors, while 30 genuine pinches remain explicit blockers.
