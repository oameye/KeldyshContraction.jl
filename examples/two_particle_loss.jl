# ## Two Body Loss
using KeldyshContraction

# ## System

# The interaction action of elastic two body scattering, is defined as
# ```math
# S_\mathrm{int} = i\Gamma \int d^d x \, [\frac{1}{2}(\bar{\phi}_+\phi_+)^2 +\frac{1}{2} (\bar{\phi}_-\phi_-)^2 -\bar{\phi}_-^2\phi_+^2]
# ```
# Above interaction can typically represent s-wave scattering of bosons.

# In the RAK basis, this gives
# ```math
# S_\mathrm{int} = i\Gamma\int d^d x \, [\frac{1}{2}(\bar{\phi}_c\bar{\phi}_q(\phi_c^2+\phi_q^2)-\phi_c\phi_q(\bar{\phi}_c^2+\bar{\phi}_q^2)+2\bar{\phi}_c\bar{\phi}_q\phi_c\phi_q)]
# ```
