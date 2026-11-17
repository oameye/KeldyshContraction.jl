"""Canonical physical identity of one spectral line for finite-width analysis."""
struct SpectralLineIdentity{S<:Statistics}
    family::FieldFamily{S}
    momentum::LinearMomentum
end

statistics(::SpectralLineIdentity{S}) where {S<:Statistics} = S
spectral_line_family(line::SpectralLineIdentity) = line.family
spectral_line_momentum(line::SpectralLineIdentity) = line.momentum

function Base.isequal(a::SpectralLineIdentity{S}, b::SpectralLineIdentity{S}) where {S}
    return isequal(a.family, b.family) && isequal(a.momentum, b.momentum)
end
Base.:(==)(a::SpectralLineIdentity, b::SpectralLineIdentity) = isequal(a, b)
function Base.hash(line::SpectralLineIdentity, h::UInt)
    return hash(line.momentum, hash(line.family, hash(SpectralLineIdentity, h)))
end
function Base.isless(a::SpectralLineIdentity{S}, b::SpectralLineIdentity{S}) where {S}
    isequal(a.family, b.family) || return isless(a.family, b.family)
    isequal(a.momentum, b.momentum) && return false
    return _linear_momentum_isless(a.momentum, b.momentum)
end

@inline function _spectral_line_identity(line::KineticLine{S}) where {S<:Statistics}
    return SpectralLineIdentity{S}(line.family, momentum(line))
end

"""
Parameters of one normalized Lorentzian/Breit--Wigner spectral line,

    A(ω) = Z Γ / ((ω - E)^2 + Γ^2/4).

`energy`, `linewidth`, and `residue` are explicit so the spectral prescription can be inspected.
The linewidth must be strictly positive. `Z = 1` is the normalized unit-residue convention.
"""
struct LorentzianSpectralData{C<:Real}
    energy::C
    linewidth::C
    residue::C

    function LorentzianSpectralData{C}(energy::C, linewidth::C, residue::C) where {C<:Real}
        linewidth > zero(C) ||
            throw(DomainError(linewidth, "spectral linewidth must be positive"))
        return new{C}(energy, linewidth, residue)
    end
end

function LorentzianSpectralData(
    energy::A, linewidth::B, residue::D
) where {A<:Real,B<:Real,D<:Real}
    C = promote_type(A, B, D)
    return LorentzianSpectralData{C}(
        convert(C, energy), convert(C, linewidth), convert(C, residue)
    )
end
function LorentzianSpectralData(energy::A, linewidth::B) where {A<:Real,B<:Real}
    C = promote_type(A, B)
    return LorentzianSpectralData{C}(convert(C, energy), convert(C, linewidth), one(C))
end

"""
    lorentzian_spectral_data_from_retarded_self_energy(ε, Σᴿ; residue=1)

Construct Lorentzian line data from an on-shell retarded self-energy using

    E = ε + Re Σᴿ,
    Γ = -2 Im Σᴿ.

`Σᴿ` is the retarded self-energy evaluated at the chosen on-shell frequency. The optional
`residue` stores an independently supplied quasiparticle residue `Z`; the default is `Z = 1`.
Frequency-dependent self-consistency and the derivative correction to `Z` are deliberately
outside this algebraic map. A nonpositive `Γ` is rejected by `LorentzianSpectralData`.
"""
function lorentzian_spectral_data_from_retarded_self_energy(
    bare_energy::A, retarded_self_energy::Complex{B}; residue=one(promote_type(A, B))
) where {A<:Real,B<:Real}
    energy = bare_energy + real(retarded_self_energy)
    linewidth = -2 * imag(retarded_self_energy)
    return LorentzianSpectralData(energy, linewidth, residue)
end

spectral_energy(data::LorentzianSpectralData) = data.energy
spectral_linewidth(data::LorentzianSpectralData) = data.linewidth
spectral_residue(data::LorentzianSpectralData) = data.residue

"""Explicit line-by-line Lorentzian spectral prescription."""
struct LorentzianSpectralModel{C<:Real,S<:Statistics}
    data::Dict{SpectralLineIdentity{S},LorentzianSpectralData{C}}
end

function LorentzianSpectralModel(
    data::AbstractDict{SpectralLineIdentity{S},LorentzianSpectralData{C}}
) where {C<:Real,S<:Statistics}
    return LorentzianSpectralModel{C,S}(Dict(data))
end

function spectral_data(
    model::LorentzianSpectralModel{C,S}, line::SpectralLineIdentity{S}
) where {C<:Real,S<:Statistics}
    haskey(model.data, line) || throw(
        ArgumentError(
            "finite-width reduction requires Lorentzian spectral data for every resolved line",
        ),
    )
    return model.data[line]
end
function spectral_energy(model::LorentzianSpectralModel, line::SpectralLineIdentity)
    return spectral_energy(spectral_data(model, line))
end
function spectral_linewidth(model::LorentzianSpectralModel, line::SpectralLineIdentity)
    return spectral_linewidth(spectral_data(model, line))
end
function spectral_residue(model::LorentzianSpectralModel, line::SpectralLineIdentity)
    return spectral_residue(spectral_data(model, line))
end

"""Evaluate the Lorentzian spectral density at frequency `ω`."""
function lorentzian_spectral_value(data::LorentzianSpectralData, ω::Real)
    Γ = spectral_linewidth(data)
    δ = ω - spectral_energy(data)
    return spectral_residue(data) * Γ / (δ^2 + Γ^2 / 4)
end

"""
Exact integrated `m`th power of one Lorentzian line,

    ∫ dω/(2π) A(ω)^m
      = binomial(2m-2,m-1) Z^m / Γ^(m-1).

For `m=1` this is the spectral normalization `Z`; for `m=2` it is `2Z²/Γ`.
"""
function lorentzian_integrated_power(data::LorentzianSpectralData, multiplicity::Integer)
    multiplicity >= 1 || throw(ArgumentError("spectral multiplicity must be positive"))
    m = Int(multiplicity)
    coefficient = binomial(2m - 2, m - 1)
    return coefficient * spectral_residue(data)^m / spectral_linewidth(data)^(m - 1)
end

spectral_normalization(data::LorentzianSpectralData) = lorentzian_integrated_power(data, 1)
spectral_squared_weight(data::LorentzianSpectralData) = lorentzian_integrated_power(data, 2)

"""One repeated spectral-line factor in an analytic Lorentzian frequency integral."""
struct LorentzianSpectralFactor{S<:Statistics}
    line::SpectralLineIdentity{S}
    multiplicity::Int

    function LorentzianSpectralFactor(
        line::SpectralLineIdentity{S}, multiplicity::Integer
    ) where {S<:Statistics}
        multiplicity >= 1 || throw(ArgumentError("spectral multiplicity must be positive"))
        return new{S}(line, Int(multiplicity))
    end
end

spectral_line(factor::LorentzianSpectralFactor) = factor.line
spectral_multiplicity(factor::LorentzianSpectralFactor) = factor.multiplicity
linewidth_exponent(factor::LorentzianSpectralFactor) = 1 - factor.multiplicity

"""
Exact analytic weight for a factorized product of Lorentzian spectral lines.

`jacobian` is `1/|det M|` for the independent affine loop-frequency transformation. Each factor
carries one physical line identity and its multiplicity. The linewidth scaling is therefore
explicit: a factor of multiplicity `m` contributes `Γ^(1-m)`.
"""
struct LorentzianSpectralWeight{S<:Statistics}
    jacobian::EnergyCoefficient
    factors::Vector{LorentzianSpectralFactor{S}}
end

function LorentzianSpectralWeight(
    jacobian::Rational, factors::AbstractVector{LorentzianSpectralFactor{S}}
) where {S<:Statistics}
    jacobian > 0 || throw(ArgumentError("spectral-frequency Jacobian must be positive"))
    canonical = collect(factors)
    sort!(canonical; by=spectral_line)
    for index in 2:length(canonical)
        spectral_line(canonical[index - 1]) == spectral_line(canonical[index]) && throw(
            ArgumentError("Lorentzian spectral factors must have unique line identities"),
        )
    end
    return LorentzianSpectralWeight{S}(convert(EnergyCoefficient, jacobian), canonical)
end

spectral_jacobian(weight::LorentzianSpectralWeight) = weight.jacobian
spectral_factors(weight::LorentzianSpectralWeight) = weight.factors

"""Structural outcome of an analytic Lorentzian finite-width reduction attempt."""
@enum LorentzianSpectralReductionKind::UInt8 begin
    LorentzianResolved = 0
    LorentzianUnsupported = 1
end

"""
Concrete result of attempting the factorized Lorentzian frequency reduction.

`LorentzianResolved` means `jacobian` and `factors` define an exact analytic spectral weight.
`LorentzianUnsupported` means the term lies outside this analytic backend; the original kinetic
term remains untouched and no numerical regulator has been introduced.
"""
struct LorentzianSpectralReduction{S<:Statistics}
    kind::LorentzianSpectralReductionKind
    jacobian::EnergyCoefficient
    factors::Vector{LorentzianSpectralFactor{S}}
end

spectral_reduction_kind(reduction::LorentzianSpectralReduction) = reduction.kind
function spectral_reduction_resolved(reduction::LorentzianSpectralReduction)
    return reduction.kind === LorentzianResolved
end

function _unsupported_lorentzian_reduction(::Type{S}) where {S<:Statistics}
    return LorentzianSpectralReduction{S}(
        LorentzianUnsupported, zero(EnergyCoefficient), LorentzianSpectralFactor{S}[]
    )
end

function _resolved_lorentzian_reduction(
    weight::LorentzianSpectralWeight{S}
) where {S<:Statistics}
    return LorentzianSpectralReduction{S}(
        LorentzianResolved, weight.jacobian, weight.factors
    )
end

function _require_resolved(reduction::LorentzianSpectralReduction)
    spectral_reduction_resolved(reduction) || throw(
        ArgumentError(
            "Lorentzian spectral reduction is unsupported for this frequency structure"
        ),
    )
    return nothing
end

function spectral_jacobian(reduction::LorentzianSpectralReduction)
    _require_resolved(reduction)
    return reduction.jacobian
end
function spectral_factors(reduction::LorentzianSpectralReduction)
    _require_resolved(reduction)
    return reduction.factors
end

"""One explicit linewidth exponent in width-aware perturbative power counting."""
struct LinewidthPower{S<:Statistics}
    line::SpectralLineIdentity{S}
    exponent::Int
end

spectral_line(power::LinewidthPower) = power.line
linewidth_exponent(power::LinewidthPower) = power.exponent

"""
Coupling order together with the linewidth powers induced by the spectral integral.

The linewidth powers are deliberately not folded into the perturbative parameter monomial. For
example, a nominal `γ²` contribution with one repeated Lorentzian line is represented as `γ²`
together with `Γ(q)^(-1)` until a separate linewidth scaling law is supplied.
"""
struct WidthAwarePowerCounting{S<:Statistics}
    parameter::ParameterMonomial
    linewidths::Vector{LinewidthPower{S}}
end

parameters(counting::WidthAwarePowerCounting) = counting.parameter
linewidth_powers(counting::WidthAwarePowerCounting) = counting.linewidths

function width_aware_power_counting(
    parameter::ParameterMonomial, weight::LorentzianSpectralWeight{S}
) where {S<:Statistics}
    linewidths = LinewidthPower{S}[]
    for factor in spectral_factors(weight)
        exponent = linewidth_exponent(factor)
        iszero(exponent) ||
            push!(linewidths, LinewidthPower(spectral_line(factor), exponent))
    end
    return WidthAwarePowerCounting{S}(parameter, linewidths)
end
function width_aware_power_counting(
    parameter::ParameterMonomial, reduction::LorentzianSpectralReduction{S}
) where {S<:Statistics}
    _require_resolved(reduction)
    linewidths = LinewidthPower{S}[]
    for factor in reduction.factors
        exponent = linewidth_exponent(factor)
        iszero(exponent) ||
            push!(linewidths, LinewidthPower(spectral_line(factor), exponent))
    end
    return WidthAwarePowerCounting{S}(parameter, linewidths)
end

"""
Evaluate an exact structural Lorentzian weight with explicit line data.

Every line in `weight` must be present in `model`; missing width data is an error rather than an
implicit sharp-shell or cutoff prescription.
"""
function evaluate_spectral_weight(
    weight::LorentzianSpectralWeight{S}, model::LorentzianSpectralModel{C,S}
) where {C<:Real,S<:Statistics}
    factors = spectral_factors(weight)
    isempty(factors) && return convert(C, spectral_jacobian(weight))

    first_factor = first(factors)
    first_value = lorentzian_integrated_power(
        spectral_data(model, spectral_line(first_factor)),
        spectral_multiplicity(first_factor),
    )
    D = typeof(first_value)
    value = convert(D, spectral_jacobian(weight)) * first_value
    for factor in Iterators.drop(factors, 1)
        contribution = lorentzian_integrated_power(
            spectral_data(model, spectral_line(factor)), spectral_multiplicity(factor)
        )
        value *= convert(D, contribution)
    end
    return value
end

function evaluate_spectral_weight(
    reduction::LorentzianSpectralReduction{S}, model::LorentzianSpectralModel{C,S}
) where {C<:Real,S<:Statistics}
    _require_resolved(reduction)
    factors = reduction.factors
    isempty(factors) && return convert(C, reduction.jacobian)

    first_factor = first(factors)
    first_value = lorentzian_integrated_power(
        spectral_data(model, spectral_line(first_factor)),
        spectral_multiplicity(first_factor),
    )
    D = typeof(first_value)
    value = convert(D, reduction.jacobian) * first_value
    for factor in Iterators.drop(factors, 1)
        contribution = lorentzian_integrated_power(
            spectral_data(model, spectral_line(factor)), spectral_multiplicity(factor)
        )
        value *= convert(D, contribution)
    end
    return value
end

"""
    lorentzian_spectral_reduction(term)

Attempt the exact analytic Lorentzian frequency reduction for an all-spectral term. Repeated
physical lines are detected from exact `(field family, routed momentum)` identity, not topology
or perturbative order.

The resolved analytic case has one independent spectral-line identity per loop frequency and a
full-rank exact frequency-routing matrix. Terms containing dispersive factors, finite Trotter
shifts, convolutions with more spectral identities than loop frequencies, or rank-deficient
routing return a concrete `LorentzianUnsupported` result. No numerical regularisation is
attempted and the original kinetic term is not modified.
"""
function lorentzian_spectral_reduction(
    term::SpectralDispersiveTerm{S}
) where {S<:Statistics}
    unsupported = _unsupported_lorentzian_reduction(S)
    lines = kinetic_lines(term.carrier)
    kinds = spectral_dispersive_kinds(term)
    all(kind -> kind === CollisionSpectral, kinds) || return unsupported
    any(line -> !iszero(regularisation_shift(line)), lines) && return unsupported

    counts = Dict{SpectralLineIdentity{S},Int}()
    for line in lines
        identity = _spectral_line_identity(line)
        counts[identity] = get(counts, identity, 0) + 1
    end
    identities = collect(keys(counts))
    sort!(identities)

    loop_indices = _loop_frequency_basis_indices(term)
    nloops = length(loop_indices)
    length(identities) == nloops || return unsupported
    if iszero(nloops)
        weight = LorentzianSpectralWeight(
            one(EnergyCoefficient), LorentzianSpectralFactor{S}[]
        )
        return _resolved_lorentzian_reduction(weight)
    end

    matrix = Matrix{EnergyCoefficient}(undef, nloops, nloops)
    for (row, identity) in pairs(identities)
        routed = spectral_line_momentum(identity)
        for (column, basis_index) in pairs(loop_indices)
            matrix[row, column] = convert(EnergyCoefficient, routed[basis_index])
        end
    end
    _exact_frequency_rank(matrix) == nloops || return unsupported
    _, determinant = _exact_inverse_and_determinant(matrix)

    factors = LorentzianSpectralFactor{S}[
        LorentzianSpectralFactor(identity, counts[identity]) for identity in identities
    ]
    weight = LorentzianSpectralWeight(abs(inv(determinant)), factors)
    return _resolved_lorentzian_reduction(weight)
end
