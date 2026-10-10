"""Physical finite sector of a quasiparticle spectral self-energy reduction."""
struct SpectralSelfEnergySector{S<:Statistics}
    parameter::ParameterMonomial
    basis::MomentumBasis
    external_momentum::MomentumVariable
    kinematic::MomentumPolynomial{ComplexRationals}
    support::FrequencySupport{S}
end

statistics(::SpectralSelfEnergySector{S}) where {S<:Statistics} = S
parameters(sector::SpectralSelfEnergySector) = sector.parameter
momentum_basis(sector::SpectralSelfEnergySector) = sector.basis
external_wigner_momentum(sector::SpectralSelfEnergySector) = sector.external_momentum
kinematic_factor(sector::SpectralSelfEnergySector) = sector.kinematic
frequency_support(sector::SpectralSelfEnergySector) = sector.support

function Base.isequal(
    a::SpectralSelfEnergySector{S}, b::SpectralSelfEnergySector{S}
) where {S<:Statistics}
    return isequal(a.parameter, b.parameter) &&
           isequal(a.basis, b.basis) &&
           isequal(a.external_momentum, b.external_momentum) &&
           isequal(a.kinematic, b.kinematic) &&
           isequal(a.support, b.support)
end
Base.:(==)(a::SpectralSelfEnergySector, b::SpectralSelfEnergySector) = isequal(a, b)
function Base.hash(sector::SpectralSelfEnergySector, h::UInt)
    h = hash(SpectralSelfEnergySector, h)
    h = hash(sector.parameter, h)
    h = hash(sector.basis, h)
    h = hash(sector.external_momentum, h)
    h = hash(sector.kinematic, h)
    return hash(sector.support, h)
end

"""
Quasiparticle spectral self-energy after exact frequency reduction, `F -> n` substitution, and
canonical dummy-loop quotient.

The regular terms represent `A_Σ = im * (Σᴿ - Σᴬ) = -2 Im Σᴿ` without applying the collision
normalization `C_n = σ I/2`. Causal-frequency blockers and unresolved finite-Trotter states are
retained explicitly, so no linewidth is manufactured when the strict quasiparticle reduction is
not defined.
"""
struct SpectralSelfEnergyKernel{
    C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext
}
    regular::Dict{SpectralSelfEnergySector{S},OccupationPolynomial{C,S}}
    blocked::Dict{TrotterFrequencyGroupKey{S},TrotterFrequencyBlockedContribution{C,S}}
    unresolved::Vector{TrotterFrequencyState{C,S,E1,E2}}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::SpectralSelfEnergyKernel{C,S,O}) where {C,S,O} = O
statistics(::SpectralSelfEnergyKernel{C,S}) where {C,S} = S
target_family(kernel::SpectralSelfEnergyKernel) = kernel.target
parameters(kernel::SpectralSelfEnergyKernel) = kernel.parameter
gradient_order(::SpectralSelfEnergyKernel{C,S,O,E1,E2,G}) where {C,S,O,E1,E2,G} = Val(G)
wigner_context(kernel::SpectralSelfEnergyKernel) = kernel.context

"""Return the finite quasiparticle occupation terms of a spectral self-energy kernel."""
spectral_self_energy_terms(kernel::SpectralSelfEnergyKernel) = kernel.regular

"""Return causal-frequency blockers retained by a spectral self-energy kernel."""
spectral_self_energy_blocked_terms(kernel::SpectralSelfEnergyKernel) = kernel.blocked

"""Return unresolved finite-Trotter states retained by a spectral self-energy kernel."""
spectral_self_energy_trotter_terms(kernel::SpectralSelfEnergyKernel) = kernel.unresolved

function Base.isempty(kernel::SpectralSelfEnergyKernel)
    return isempty(kernel.regular) && isempty(kernel.blocked) && isempty(kernel.unresolved)
end

function _spectral_self_energy_sector(sector)
    S = statistics(sector)
    return SpectralSelfEnergySector{S}(
        parameters(sector),
        momentum_basis(sector),
        external_wigner_momentum(sector),
        kinematic_factor(sector),
        frequency_support(sector),
    )
end

function _spectral_self_energy_frequency_problem(
    Σ::KineticSelfEnergy{C,S,O,E1,E2,G,Ctx}
) where {C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    spectral = _spectral_dispersive_expression(spectral_self_energy(Σ))
    empty = SpectralDispersiveExpression{C,S,E1,E2,G,Ctx}(wigner_context(Σ))
    return SpectralDispersiveCollision{C,S,O,E1,E2,G,Ctx}(
        spectral, empty, target_family(Σ), parameters(Σ), wigner_context(Σ)
    )
end

"""
    spectral_self_energy_kernel(Σ)

Reduce the generated spectral self-energy `A_Σ = im * (Σᴿ - Σᴬ)` to the strict-quasiparticle
occupation representation while preserving exact canonical loop routing, shell/PV support,
parameter provenance, and unresolved causal/Trotter states.

This reuses the collision frequency reducer and physical dummy-loop quotient only as exact algebra
engines. The statistical lowering is `occupation_substitute`, not
`occupation_collision_polynomial`, so a self-energy linewidth never acquires the collision
normalization `σ/2`.
"""
function spectral_self_energy_kernel(
    Σ::KineticSelfEnergy{C,S,O,E1,E2,G,Ctx}
) where {C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    reduced = reduce_frequency_collision(_spectral_self_energy_frequency_problem(Σ))
    occupation = _occupation_reduced_expression(reduced, occupation_substitute)
    quotient = quotient_loop_momenta(occupation)

    D = promote_type(C, ComplexRationals, Rational{Int})
    regular = Dict{SpectralSelfEnergySector{S},OccupationPolynomial{D,S}}()
    for (sector, polynomial) in loop_quotient_terms(quotient)
        regular[_spectral_self_energy_sector(sector)] = polynomial
    end

    return SpectralSelfEnergyKernel{D,S,O,E1,E2,G,Ctx}(
        regular,
        copy(reduced_blocked_terms(reduced)),
        copy(reduced_trotter_terms(reduced)),
        target_family(Σ),
        parameters(Σ),
        wigner_context(Σ),
    )
end
