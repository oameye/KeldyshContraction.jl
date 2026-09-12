# ## Fermionic p-wave Two-Body Loss
using KeldyshContraction
using KeldyshContraction: Regularisation

# ## Contour interaction

# For identical spinless fermions, local two-body loss is derivative coupled. For one
# Cartesian channel take
# ```math
# J_x = \psi\,\partial_x\psi.
# ```
# The asymmetric fermionic Larkin-Ovchinnikov rotation is built into the branch fields
# `One` and `Two`. We construct the finite-Trotter contour dissipator directly so that the
# same equal-time provenance used by the bosonic loss problem survives the kinetic reduction.

@qfields ψ::Fermion
ψ₁ = ψ[One]
ψ₂ = ψ[Two]
bψ₁ = bar(ψ₁)
bψ₂ = bar(ψ₂)
plus = Regularisation.Plus
minus = Regularisation.Minus

ψplus_minus = ψ₁(minus) + ψ₂(minus)
∂ψplus_minus = partial(ψ₁(minus), :x) + partial(ψ₂(minus), :x)
ψplus_plus = ψ₁(plus) + ψ₂(plus)
∂ψplus_plus = partial(ψ₁(plus), :x) + partial(ψ₂(plus), :x)
ψminus_plus = ψ₁(plus) - ψ₂(plus)
∂ψminus_plus = partial(ψ₁(plus), :x) - partial(ψ₂(plus), :x)

bψplus = bψ₁ + bψ₂
bψminus = bψ₂ - bψ₁
∂bψplus = partial(bψ₁, :x) + partial(bψ₂, :x)
∂bψminus = partial(bψ₂, :x) - partial(bψ₁, :x)

Jplus_minus = ψplus_minus * ∂ψplus_minus
Jplus_plus = ψplus_plus * ∂ψplus_plus
Jminus_plus = ψminus_plus * ∂ψminus_plus
Jplus_dagger = ∂bψplus * bψplus
Jminus_dagger = ∂bψminus * bψminus

loss =
    (1 // 8) *
    im *
    (
        (Jplus_dagger - Jminus_dagger) * Jplus_minus -
        (Jplus_plus - Jminus_plus) * Jminus_dagger
    )

L_int = InteractionLagrangian(loss, :γp)

# ## End-to-end kinetic kernel

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

# Antisymmetry combines the ordered derivative pieces into the relative p-wave matrix
# element. The resulting first-order occupation kernel is
# ```math
# C_n^{\gamma_p}(k)
# = -\gamma_p\int_q (k_x-q_x)^2 n_k n_q.
# ```
# No final-state Pauli-blocking factor appears because this channel removes a pair; Fermi
# statistics enters through the antisymmetric derivative matrix element.
kernel
