# ## Elastic Two Body Scattering
using KeldyshContraction

# ## System

# The interaction action of elastic two body scattering, is defined as
# ```math
# S_\mathrm{int} = -\frac{g}{2} \int d^d x \, [(\bar{\phi}_+\phi_+)^2 - (\bar{\phi}_-\phi_-)^2 ]
# ```
# Above interaction can typically represent s-wave scattering of bosons. Here, the overall
# minus sign comes from the relation of the Lagrangian to the Hamiltonian.

# In the RAK basis, this gives
# ```math
# S_\mathrm{int} = -\frac{g}{2} \int d^d x \, [(\bar{\phi}_c\bar{\phi}_q\phi_c\phi_c)
# +(\bar{\phi}_c\bar{\phi}_q\phi_q\phi_q) + c.c.]
# ```

# Let us represent this quantum and classical field inside the `KeldyshContraction` package
# and define the interaction Lagrangian:

@qfields c::Destroy(Classical) q::Destroy(Quantum)
elasctic2boson = -(0.5 * (c^2 + q^2) * c' * q' + 0.5 * c * q * ((c')^2 + (q')^2))
L_int = InteractionLagrangian(elasctic2boson)

# A good check if the interaction Lagrangian is a valid physical process, is to check if the
# normalization identity $Z=1$ holds. We can do this perturbatively in $g$ by expanding
# $\exp(i S_\mathrm{int})$  and showing the average of the linear part of the system is zero
# ```math
# \langle S_\mathrm{int}\rangle =  \langle S_\mathrm{int}^2\rangle  =\ldots = 0
# ```
# In computing the average, one performs Wick contractions to describe the average in terms
# of the two-point correlators of the linear part of the system.

# In the package we can do this as follows:
KeldyshContraction._wick_contraction(elasctic2boson; simplify=false)

# However, to show that these diagrams cancel out, we need to apply to condition $G^R = - G^A$.
# Inside the package we do this by
KeldyshContraction._wick_contraction(elasctic2boson; simplify=true)

# Similarly, we can compute the next orders.

# ## First order Green's function

# To compute the two point Green's functions of the interacting system, we can apply
#  self-consistent perturbation theory. By again expanding in $g$ we can write
# ```math
# \begin{aligned}
# & i G^{\mu \nu}\left(x_1, x_2\right)=\int \mathcal{D}\left[\phi_c, \bar{\phi}_c, \phi_q, \bar{\phi}_q\right] \phi_\mu\left(x_1\right) \bar{\phi}_\nu\left(x_2\right) \sum_{k=0}^{\infty} \frac{i^k S_{\mathrm{int}}^k}{k!} e^{i S_0}= \\
# & \quad=i G_0^{\mu \nu}\left(x_1, x_2\right)+i \int d^d y d t_y\left\langle\phi_\mu\left(x_1\right) \bar{\phi}_\nu\left(x_2\right) \mathcal{L}_{\mathrm{int}}(y)\right\rangle_0+\sum_{k=2}^{\infty}\left\langle\phi_\mu\left(x_1\right) \bar{\phi}_\nu\left(x_2\right) \frac{i^k S_{\mathrm{int}}^k}{k!}\right\rangle_0 .
# \end{aligned}
# ```

# So we can compute the first order Green's function correction $G_{(1)}$ by computing
# the Wick contraction of the interaction Lagrangian

GF = DressedPropagator(L_int)

# Here, the simplification of the advanced to retarded propagator is done internally.

# ## Self-Energy

# Often we are interested in the self-energy of the system, which is defined as
# the set of irreducible diagrams. Inside the package we can compute these to an order $g$ by

Σ = SelfEnergy(GF)

# ## Transport

# The self-energy can be used to compute derive a kinetic equation for the system.
# In doing this one compute the so-called collision integral, which is given by
# ```math
# I _\mathrm{coll}= i Σ^K(x, p) +  F (x, p) (Σ^R(x, p)-Σ^A(x, p)).
# ```
# Here, $F$ is the bosonic distribution function of the system.
# However, from above calculation we find that ``i Σ^K(x, p) = 0`` and ``Σ^R(x, p)=Σ^A(x, p)``,
# such that the collision integral has not contribution at first order.

# ## Second order

# In second order, we have many additional terms for the dressed propagator.
# These involve now 5 propagators:

GF = DressedPropagator(L_int, 2)

# However, not all of them contribute to the second order self-energy. Indeed, many terms
# (diagrams) involve only first order self-energy corrections and are thus reducible.
# Instead, we need to separate the reducible and irreducible diagrams. We can separate them
# by looking at the multiplicity of the edges in the diagrams.

topology_dict = topologies(GF.keldysh)

# This gives us three distinct topologies, which we can identify by the multiplicity of the edges.

topology_dict[[2]]

# The topology involving only one edge is the reducible diagram, which will not contribute to the self-energy in second order. Indeed, internally we only consider the irreducible diagrams.

Σ = SelfEnergy(GF, 2)
