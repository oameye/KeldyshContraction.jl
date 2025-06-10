# ## Two Body Loss

using KeldyshContraction
using KeldyshContraction: Classical, Quantum, Plus, Minus

@qfields c::Destroy(Classical) q::Destroy(Quantum)
@syms Γ::Real

loss2boson =
    im*(
        0.5 * c' * q' * (c(Minus) * c(Minus) + q(Minus) * q(Minus)) -
        0.5 * c(Plus) * q(Plus) * (c' * c' + q' * q') #+
        # c' * q' * (c(Plus) * q(Plus) + c(Minus) * q(Minus))
    )
L_int = InteractionLagrangian(loss2boson)

GF = DressedPropagator(L_int)

Σ = SelfEnergy(GF)

using KeldyshContraction: Diagram, Diagrams, Bulk, Edge, In, Out
import KeldyshContraction as KC
using Combinatorics

order = 2
in_out = c(Out()) * q'(In())
l = length(L_int.lagrangian)

E = KC.number_of_propagators(L_int) * order + 1 # +1 for in_out
dict = OrderedDict()
keys = []
prefactor = -1 * im * im^order / factorial(order)

for coefficients in Combinatorics.multiexponents(l, order)
    idxs = KC.expand_coefficients(coefficients) # will be of length order
    mult = Combinatorics.multinomial(coefficients...)
    qmul = mult * prod(L_int(i).lagrangian.arguments[j] for (i, j) in pairs(idxs))

    term = prefactor * in_out * qmul
    diagrams = wick_contraction(term; simplify=true)
    KC._simplify_prefactors!(diagrams)
    dict[term] = diagrams
    push!(keys, term)
end
dict

GF = DressedPropagator(L_int; order=2)

Σ = SelfEnergy(GF; order=2)
