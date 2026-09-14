# # Fermionic p-wave scattering and two-body loss
#
# Identical spinless fermions require an odd-wave two-particle operator. For one Cartesian
# channel take `P_x = ψ ∂xψ`. We include both the coherent interaction
# `H_int = (g_p/2)P_x†P_x` and the jump operator `L_x = √γ_p P_x`.

using KeldyshContraction
using KeldyshContraction: Regularisation

@qfields ψ::Fermion
ψ1, ψ2 = ψ[One], ψ[Two]
bψ1, bψ2 = bar(ψ1), bar(ψ2)

# ## Coherent p-wave interaction

ψplus = ψ1 + ψ2
ψminus = ψ1 - ψ2
∂ψplus = partial(ψ1, :x) + partial(ψ2, :x)
∂ψminus = partial(ψ1, :x) - partial(ψ2, :x)
bψplus = bψ1 + bψ2
bψminus = bψ2 - bψ1
∂bψplus = partial(bψ1, :x) + partial(bψ2, :x)
∂bψminus = partial(bψ2, :x) - partial(bψ1, :x)

Pplus = ψplus * ∂ψplus
Pminus = ψminus * ∂ψminus
Pplus_dagger = ∂bψplus * bψplus
Pminus_dagger = ∂bψminus * bψminus

elastic = -(1 // 8) * (Pplus_dagger * Pplus - Pminus_dagger * Pminus)

# ## Finite-Trotter loss vertex

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
L = InteractionLagrangian(elastic, gp) + InteractionLagrangian(loss, γp)

# ## Self-energy to collision kernel

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

G1 = DressedPropagator(L, Val(1), Val(3); preserve_regularisation=true)
Cγp = kinetic_kernel(G1[γp])

# Antisymmetry combines the derivative vertices into the relative p-wave matrix element,
# ```math
# C_n^{(\gamma_p)}(k)
# =-\gamma_p\int_q (k_x-q_x)^2 n_kn_q.
# ```
Cγp

G2 = DressedPropagator(L, Val(2), Val(5); preserve_regularisation=true)
Cgp2 = kinetic_kernel(G2[gp^2])

# The elastic sector has the fermionic gain--loss structure
# ```math
# \frac{g_p^2}{8}(k_x-p_x)^2(q_x-r_x)^2
# \big[(1-n_k)(1-n_p)n_qn_r
#      -n_kn_p(1-n_q)(1-n_r)\big]\,\delta(\Delta E).
# ```
# The same shared compiler also retains the mixed `g_p γ_p` principal-value sector and leaves
# genuine strict-quasiparticle pinches explicit rather than assigning them a finite value.
Cgp2
