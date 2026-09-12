const EnergyCoefficient = MomentumCoefficient

"""One symbolic quasiparticle energy `ε_family(momentum)` in an exact momentum basis."""
struct DispersionAtom{S<:Statistics}
    family::FieldFamily{S}
    momentum::LinearMomentum
end

statistics(::DispersionAtom{S}) where {S} = S
momentum(atom::DispersionAtom) = atom.momentum

function Base.isequal(a::DispersionAtom{S}, b::DispersionAtom{S}) where {S<:Statistics}
    return isequal(a.family, b.family) && isequal(a.momentum, b.momentum)
end
Base.:(==)(a::DispersionAtom, b::DispersionAtom) = isequal(a, b)
function Base.hash(atom::DispersionAtom, h::UInt)
    return hash(atom.momentum, hash(atom.family, hash(DispersionAtom, h)))
end
function Base.isless(a::DispersionAtom{S}, b::DispersionAtom{S}) where {S<:Statistics}
    isequal(a.family, b.family) || return isless(a.family, b.family)
    isequal(a.momentum, b.momentum) && return false
    return _linear_momentum_isless(a.momentum, b.momentum)
end

"""
Canonical exact rational linear combination of symbolic dispersion atoms.

`basis_size` is retained even for the zero form so dimension checks do not depend on runtime
term values.
"""
struct EnergyForm{S<:Statistics}
    terms::Vector{Pair{DispersionAtom{S},EnergyCoefficient}}
    basis_size::Int16

    function EnergyForm{S}(
        terms::Vector{Pair{DispersionAtom{S},EnergyCoefficient}},
        basis_size::Int16,
        ::Val{:raw},
    ) where {S<:Statistics}
        return new{S}(terms, basis_size)
    end
end

function EnergyForm{S}(basis_size::Integer) where {S<:Statistics}
    0 <= basis_size <= typemax(Int16) || throw(ArgumentError("invalid energy basis size"))
    return EnergyForm{S}(
        Pair{DispersionAtom{S},EnergyCoefficient}[], convert(Int16, basis_size), Val(:raw)
    )
end

function _canonical_energy_form(
    terms::Vector{Pair{DispersionAtom{S},EnergyCoefficient}}, basis_size::Int16
) where {S<:Statistics}
    for (atom, _) in terms
        length(momentum(atom)) == basis_size ||
            throw(DimensionMismatch("dispersion atom uses a different momentum basis size"))
    end
    sort!(terms; by=first)
    out = Pair{DispersionAtom{S},EnergyCoefficient}[]
    sizehint!(out, length(terms))
    for (atom, coefficient) in terms
        iszero(coefficient) && continue
        if !isempty(out) && isequal(first(out[end]), atom)
            combined = last(out[end]) + coefficient
            pop!(out)
            iszero(combined) || push!(out, atom => combined)
        else
            push!(out, atom => coefficient)
        end
    end
    return EnergyForm{S}(out, basis_size, Val(:raw))
end

function EnergyForm(
    basis_size::Integer, terms::AbstractVector{<:Pair{DispersionAtom{S},<:Rational}}
) where {S<:Statistics}
    0 <= basis_size <= typemax(Int16) || throw(ArgumentError("invalid energy basis size"))
    owned = Pair{DispersionAtom{S},EnergyCoefficient}[
        atom => convert(EnergyCoefficient, coefficient) for (atom, coefficient) in terms
    ]
    return _canonical_energy_form(owned, convert(Int16, basis_size))
end

function EnergyForm(
    atom::DispersionAtom{S}, coefficient::Rational=1 // 1
) where {S<:Statistics}
    return EnergyForm(length(momentum(atom)), [atom => coefficient])
end

energy_terms(form::EnergyForm) = form.terms
energy_basis_size(form::EnergyForm) = Int(form.basis_size)
Base.length(form::EnergyForm) = length(form.terms)
Base.isempty(form::EnergyForm) = isempty(form.terms)
Base.iszero(form::EnergyForm) = isempty(form.terms)
Base.iterate(form::EnergyForm) = iterate(form.terms)
Base.iterate(form::EnergyForm, state) = iterate(form.terms, state)
Base.eltype(::Type{EnergyForm{S}}) where {S} = Pair{DispersionAtom{S},EnergyCoefficient}
Base.zero(form::EnergyForm{S}) where {S} = EnergyForm{S}(energy_basis_size(form))

function Base.isequal(a::EnergyForm{S}, b::EnergyForm{S}) where {S<:Statistics}
    return a.basis_size == b.basis_size && isequal(a.terms, b.terms)
end
Base.:(==)(a::EnergyForm, b::EnergyForm) = isequal(a, b)
function Base.hash(form::EnergyForm, h::UInt)
    return hash(form.terms, hash(form.basis_size, hash(EnergyForm, h)))
end

function Base.isless(a::EnergyForm{S}, b::EnergyForm{S}) where {S<:Statistics}
    a.basis_size == b.basis_size || return a.basis_size < b.basis_size
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        a_atom, a_coefficient = a.terms[i]
        b_atom, b_coefficient = b.terms[i]
        isequal(a_atom, b_atom) || return isless(a_atom, b_atom)
        a_coefficient == b_coefficient || return a_coefficient < b_coefficient
    end
    return length(a) < length(b)
end

function _check_energy_dimensions(a::EnergyForm, b::EnergyForm)
    a.basis_size == b.basis_size || throw(DimensionMismatch("energy basis sizes differ"))
    return nothing
end

function Base.:+(a::EnergyForm{S}, b::EnergyForm{S}) where {S<:Statistics}
    _check_energy_dimensions(a, b)
    terms = Pair{DispersionAtom{S},EnergyCoefficient}[]
    sizehint!(terms, length(a) + length(b))
    append!(terms, a.terms)
    append!(terms, b.terms)
    return _canonical_energy_form(terms, a.basis_size)
end
Base.:-(form::EnergyForm) = (-1 // 1) * form
Base.:-(a::EnergyForm, b::EnergyForm) = a + (-b)

function Base.:*(coefficient::Rational, form::EnergyForm{S}) where {S<:Statistics}
    c = convert(EnergyCoefficient, coefficient)
    terms = Pair{DispersionAtom{S},EnergyCoefficient}[
        atom => c * value for (atom, value) in form
    ]
    return _canonical_energy_form(terms, form.basis_size)
end
Base.:*(coefficient::Integer, form::EnergyForm) = (coefficient // 1) * form
Base.:*(form::EnergyForm, coefficient::Union{Integer,Rational}) = coefficient * form

function _canonical_energy_scale(form::EnergyForm)
    iszero(form) && throw(ArgumentError("zero energy form has no support scale"))
    scale = last(first(form.terms))
    return inv(scale) * form, scale
end

"""Canonical energy-shell support denominator with leading exact coefficient one."""
struct EnergyShell{S<:Statistics}
    energy::EnergyForm{S}

    function EnergyShell{S}(energy::EnergyForm{S}, ::Val{:canonical}) where {S}
        return new{S}(energy)
    end
end

"""
Return `(support, factor)` for `δ(E)` in canonical energy orientation and scale.

If `E = c E₀` with the first canonical coefficient of `E₀` equal to one, the returned factor
is `1/abs(c)` because `δ(c E₀) = δ(E₀)/abs(c)`.
"""
function energy_shell(form::EnergyForm{S}) where {S<:Statistics}
    canonical, scale = _canonical_energy_scale(form)
    return EnergyShell{S}(canonical, Val(:canonical)), inv(abs(scale))
end

Base.isequal(a::EnergyShell{S}, b::EnergyShell{S}) where {S} = isequal(a.energy, b.energy)
Base.:(==)(a::EnergyShell, b::EnergyShell) = isequal(a, b)
Base.hash(shell::EnergyShell, h::UInt) = hash(shell.energy, hash(EnergyShell, h))
Base.isless(a::EnergyShell, b::EnergyShell) = isless(a.energy, b.energy)

"""Canonical denominator stored by a principal-value support factor."""
struct PrincipalValueSupport{S<:Statistics}
    energy::EnergyForm{S}

    function PrincipalValueSupport{S}(energy::EnergyForm{S}, ::Val{:canonical}) where {S}
        return new{S}(energy)
    end
end

"""
Return `(support, factor)` for `PV(1/E)` in canonical denominator orientation and scale.

If `E = c E₀`, the returned factor is `1/c`; unlike a delta function, reversing the denominator
therefore reverses the principal-value contribution's sign.
"""
function principal_value_support(form::EnergyForm{S}) where {S<:Statistics}
    canonical, scale = _canonical_energy_scale(form)
    return PrincipalValueSupport{S}(canonical, Val(:canonical)), inv(scale)
end

function Base.isequal(a::PrincipalValueSupport{S}, b::PrincipalValueSupport{S}) where {S}
    return isequal(a.energy, b.energy)
end
Base.:(==)(a::PrincipalValueSupport, b::PrincipalValueSupport) = isequal(a, b)
function Base.hash(support::PrincipalValueSupport, h::UInt)
    return hash(support.energy, hash(PrincipalValueSupport, h))
end
Base.isless(a::PrincipalValueSupport, b::PrincipalValueSupport) = isless(a.energy, b.energy)

"""Canonical collection of residual shell and principal-value frequency supports."""
struct FrequencySupport{S<:Statistics}
    shells::Vector{EnergyShell{S}}
    principal_values::Vector{PrincipalValueSupport{S}}

    function FrequencySupport{S}(
        shells::Vector{EnergyShell{S}},
        principal_values::Vector{PrincipalValueSupport{S}},
        ::Val{:canonical},
    ) where {S<:Statistics}
        return new{S}(shells, principal_values)
    end
end

function FrequencySupport(
    shells::AbstractVector{EnergyShell{S}},
    principal_values::AbstractVector{PrincipalValueSupport{S}},
) where {S<:Statistics}
    canonical_shells = sort!(collect(shells))
    canonical_pv = sort!(collect(principal_values))
    return FrequencySupport{S}(canonical_shells, canonical_pv, Val(:canonical))
end

function Base.isequal(a::FrequencySupport{S}, b::FrequencySupport{S}) where {S}
    return isequal(a.shells, b.shells) && isequal(a.principal_values, b.principal_values)
end
Base.:(==)(a::FrequencySupport, b::FrequencySupport) = isequal(a, b)
function Base.hash(support::FrequencySupport, h::UInt)
    return hash(support.principal_values, hash(support.shells, hash(FrequencySupport, h)))
end
