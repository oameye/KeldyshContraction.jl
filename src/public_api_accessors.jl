@doc """
    field_families(x)

Return the physical field families carried by an interaction.
""" field_families

@doc """
    target_family(x)

Return the physical field family selected for the external two-point function or kinetic
result.
""" target_family

@doc """
    parameters(x)

Return the perturbative parameter monomial, or the available monomials for a multi-process
result.
""" parameters

@doc """
    matrix(x)

Assemble the statistics-specific Keldysh matrix from the semantic retarded, advanced, and
Keldysh components of a propagator or self-energy.
""" matrix

@doc """
    order(x)

Return the perturbative order represented by `x`.
""" order

@doc """
    statistics(x)

Return the field statistics represented by `x`.
""" statistics

@doc """
    gradient_order(x)

Return the retained Wigner-gradient order as a `Val`. The current kinetic compiler supports
only the homogeneous zeroth-gradient collision problem.
""" gradient_order

@doc """
    wigner_context(x)

Return the Wigner context carried by a Wigner-space or kinetic result.
""" wigner_context

@doc """
    topologies(diagrams)

Group diagrams by their canonical uncolored topology signature. This is an advanced inspection
helper; physical diagram identity still retains field, Keldysh, regularisation, and propagator
information beyond the topology label.
""" topologies

@doc """
    wigner_transform(x; gradient_order=Val(0))

Transform a Fourier-space propagator or self-energy to homogeneous Wigner variables while
preserving exact momentum routing, derivative kinematics, target provenance, and finite
equal-time regularisation. The current kinetic compiler supports only `gradient_order=Val(0)`.
""" wigner_transform

"""
    keldysh_component(x)

Return the Keldysh component of a propagator, self-energy, Fourier/Wigner result, or kinetic
self-energy returned by the public workflow.
"""
keldysh_component(x::DressedPropagator) = x.keldysh
keldysh_component(x::SelfEnergy) = x.keldysh
keldysh_component(x::FourierDressedPropagator) = x.keldysh
keldysh_component(x::FourierSelfEnergy) = x.keldysh
keldysh_component(x::WignerDressedPropagator) = x.keldysh
keldysh_component(x::WignerSelfEnergy) = x.keldysh
keldysh_component(x::KineticSelfEnergy) = x.keldysh

"""
    retarded_component(x)

Return the retarded component of a propagator, self-energy, Fourier/Wigner result, or kinetic
self-energy returned by the public workflow.
"""
retarded_component(x::DressedPropagator) = x.retarded
retarded_component(x::SelfEnergy) = x.retarded
retarded_component(x::FourierDressedPropagator) = x.retarded
retarded_component(x::FourierSelfEnergy) = x.retarded
retarded_component(x::WignerDressedPropagator) = x.retarded
retarded_component(x::WignerSelfEnergy) = x.retarded
retarded_component(x::KineticSelfEnergy) = x.retarded

"""
    advanced_component(x)

Return the advanced component of a propagator, self-energy, Fourier/Wigner result, or kinetic
self-energy returned by the public workflow.
"""
advanced_component(x::DressedPropagator) = x.advanced
advanced_component(x::SelfEnergy) = x.advanced
advanced_component(x::FourierDressedPropagator) = x.advanced
advanced_component(x::FourierSelfEnergy) = x.advanced
advanced_component(x::WignerDressedPropagator) = x.advanced
advanced_component(x::WignerSelfEnergy) = x.advanced
advanced_component(x::KineticSelfEnergy) = x.advanced

@doc """
    collision_offset(collision)

Return the distribution-independent part of a collision expression. For
`OffShellCollisionExpression` this is the exact full-frequency contribution
`I₀ = iΣᴷ`; no quasiparticle shell, occupation substitution, or finite-width prescription has
been applied.
""" collision_offset

@doc """
    collision_distribution_coefficient(collision)

Return the coefficient of the external statistical distribution in a collision expression. For
`OffShellCollisionExpression` this is the exact full-frequency contribution `I₁ = -A_Σ`, so
`I_coll = I₀ + F_target(k) I₁`. The external distribution is not absorbed into the internal-line
IR and the quasiparticle reduction remains a downstream operation.
""" collision_distribution_coefficient

@doc """
    kinematic_factor(sector)

Return the exact momentum polynomial multiplying a reduced collision sector.
""" kinematic_factor

@doc """
    frequency_support(sector)

Return the shell and principal-value frequency support carried by a reduced collision sector.
""" frequency_support

@doc """
    reduced_regular_terms(result)

Return the finite shell and principal-value terms after causal frequency reduction.
""" reduced_regular_terms

@doc """
    reduced_blocked_terms(result)

Return causal frequency structures that do not have a finite value in the strict reduction,
including genuine pinch singularities.
""" reduced_blocked_terms

@doc """
    reduced_trotter_terms(result)

Return finite-Trotter states that remain unresolved after structural equal-time reduction.
""" reduced_trotter_terms

@doc """
    collision_kernel_terms(kernel)

Return the canonical sector-to-occupation-polynomial mapping of a final `CollisionKernel`.
""" collision_kernel_terms
