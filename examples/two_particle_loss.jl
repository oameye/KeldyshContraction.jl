# ## Two Body Loss

using KeldyshContraction
using KeldyshContraction: Classical, Quantum, Plus, Minus

@qfields ϕ::Destroy(Classical) ψ::Destroy(Quantum)
@syms Γ::Real

loss2boson =
    im*(
        0.5 * ϕ' * ψ' * (ϕ(Minus) * ϕ(Minus) + ψ(Minus) * ψ(Minus)) -
        0.5 * ϕ(Plus) * ψ(Plus) * (ϕ' * ϕ' + ψ' * ψ') +
        ϕ' * ψ' * (ϕ(Plus) * ψ(Plus) + ϕ(Minus) * ψ(Minus))
    )
L_int = InteractionLagrangian(loss2boson)

GF = DressedPropagator(L_int)

Σ = SelfEnergy(GF)

GF = DressedPropagator(L_int; order=2)

Σ = SelfEnergy(GF; order=2)
