
using KeldyshContraction

@qfields c::Destroy(Classical) q::Destroy(Quantum)
elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
L_int = InteractionLagrangian(elasctic2boson)

GF = DressedPropagator(L_int, 2)

Σ = SelfEnergy(GF, 2)

Σk = wigner_transform(Σ)
diagram = first(first(Σk.retarded))
KeldyshContraction.imaginary_part(diagram)

dict = KeldyshContraction.imaginary_part(Σk.retarded)
repr(dict[[3]])
