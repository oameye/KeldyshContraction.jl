"""
Spectral self-energy kernel after evaluating only its occupation factors on a supplied background.

The routed momentum basis, kinematic factor, shell/PV support, perturbative parameter, and any
causal/Trotter provenance remain unchanged. The stored scalar is therefore the coefficient of the
still-explicit routed momentum integral; this type does not perform momentum integration or turn a
linewidth integrand into a scalar `Γ`.
"""
struct BackgroundSpectralSelfEnergyKernel{
    V<:Number,C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext
}
    regular::Dict{SpectralSelfEnergySector{S},V}
    blocked::Dict{TrotterFrequencyGroupKey{S},TrotterFrequencyBlockedContribution{C,S}}
    unresolved::Vector{TrotterFrequencyState{C,S,E1,E2}}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::BackgroundSpectralSelfEnergyKernel{V,C,S,O}) where {V,C,S,O} = O
statistics(::BackgroundSpectralSelfEnergyKernel{V,C,S}) where {V,C,S} = S
target_family(kernel::BackgroundSpectralSelfEnergyKernel) = kernel.target
parameters(kernel::BackgroundSpectralSelfEnergyKernel) = kernel.parameter
function gradient_order(
    ::BackgroundSpectralSelfEnergyKernel{V,C,S,O,E1,E2,G}
) where {V,C,S,O,E1,E2,G}
    return Val(G)
end
wigner_context(kernel::BackgroundSpectralSelfEnergyKernel) = kernel.context

"""Return background-evaluated regular spectral-self-energy integrand coefficients."""
function background_spectral_self_energy_terms(kernel::BackgroundSpectralSelfEnergyKernel)
    return kernel.regular
end

function spectral_self_energy_blocked_terms(kernel::BackgroundSpectralSelfEnergyKernel)
    return kernel.blocked
end
function spectral_self_energy_trotter_terms(kernel::BackgroundSpectralSelfEnergyKernel)
    return kernel.unresolved
end

function Base.isempty(kernel::BackgroundSpectralSelfEnergyKernel)
    return isempty(kernel.regular) && isempty(kernel.blocked) && isempty(kernel.unresolved)
end

"""
    evaluate_spectral_self_energy_background(kernel, background)

Evaluate the exact occupation polynomial of every regular spectral-self-energy sector on a typed
`OccupationBackground`. Zero regular sectors are removed canonically. Momentum routing,
derivative kinematics, shell/PV support, and unresolved causal/Trotter branches are preserved.

This operation deliberately stops before loop-momentum integration. For example, evaluating the
first-order two-body-loss kernel `4 n_q` on a background with `n_q = m` yields the routed integrand
coefficient `4m`; obtaining `Γ = 4γ∫_q n_q` still requires an explicit user-supplied momentum
integration/closure.
"""
function evaluate_spectral_self_energy_background(
    kernel::SpectralSelfEnergyKernel{C,S,O,E1,E2,G,Ctx}, background::OccupationBackground{V}
) where {C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext,V<:Number}
    D = promote_type(C, V)
    regular = Dict{SpectralSelfEnergySector{S},D}()
    for (sector, polynomial) in spectral_self_energy_terms(kernel)
        value = convert(D, evaluate_occupation_polynomial(polynomial, background))
        iszero(value) || (regular[sector] = value)
    end
    return BackgroundSpectralSelfEnergyKernel{D,C,S,O,E1,E2,G,Ctx}(
        regular,
        copy(spectral_self_energy_blocked_terms(kernel)),
        copy(spectral_self_energy_trotter_terms(kernel)),
        target_family(kernel),
        parameters(kernel),
        wigner_context(kernel),
    )
end
