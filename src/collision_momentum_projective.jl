# Projective normalization of routed derivative factors at the physical loop quotient.
#
# A MomentumComponent stores an exact linear form such as q or -q. Scalar multiples
# of the same linear form are algebraically the same polynomial generator; the scalar
# belongs in the MomentumPolynomial coefficient. Canonicalizing that scalar before a
# dummy-loop gauge is chosen prevents coordinate orientation from erasing derivative
# signs.

function _projective_kinematic_monomial(monomial::MomentumMonomial)
    components = MomentumComponent[]
    sizehint!(components, length(monomial))
    factor = one(MomentumCoefficient)

    for component in monomial
        pivot = zero(MomentumCoefficient)
        for coefficient in component.momentum.coefficients
            if !iszero(coefficient)
                pivot = coefficient
                break
            end
        end
        iszero(pivot) && return MomentumMonomial(), zero(MomentumCoefficient), false

        normalized_momentum = inv(pivot) * component.momentum
        push!(components, MomentumComponent(normalized_momentum, component.axis))
        factor *= pivot
    end

    return MomentumMonomial(components), factor, true
end

# Derivative lowering uses this exact coefficient domain throughout the kinetic path.
# Put each routed linear factor in projective normal form at construction time so the
# later loop quotient sees q_x and -q_x as one generator with coefficients +1 and -1,
# rather than as two coordinate spellings whose signs can be absorbed by different
# dummy-loop gauges.
function MomentumPolynomial(
    monomial::MomentumMonomial, coefficient::ComplexRationals
)::MomentumPolynomial{ComplexRationals}
    normalized, factor, nonzero = _projective_kinematic_monomial(monomial)
    (nonzero && !iszero(coefficient)) || return MomentumPolynomial{ComplexRationals}()
    normalized_coefficient = coefficient * convert(ComplexRationals, factor)
    return MomentumPolynomial{ComplexRationals}([normalized => normalized_coefficient])
end

function transform_loop_momenta(
    polynomial::MomentumPolynomial{ComplexRationals}, transform::LoopMomentumTransform
)
    terms = Pair{MomentumMonomial,ComplexRationals}[]
    sizehint!(terms, length(polynomial))

    for (monomial, coefficient) in polynomial
        transformed = transform_loop_momenta(monomial, transform)
        normalized, factor, nonzero = _projective_kinematic_monomial(transformed)
        nonzero || continue
        push!(terms, normalized => coefficient * convert(ComplexRationals, factor))
    end

    return MomentumPolynomial{ComplexRationals}(terms)
end

# Signed dummy-loop reflections are part of the quotient group. If one such
# reflection leaves the occupation monomial and frequency support unchanged while
# multiplying the derivative kinematics by c != 1, the complete atom satisfies
# I = c I and therefore vanishes exactly. This does not assume n(q) = n(-q): atoms
# containing the reflected occupation momentum fail the invariance test and survive.
function _vanishes_under_loop_reflection(
    sector::CollisionKernelSector{S}, monomial::OccupationMonomial{S}
) where {S<:Statistics}
    basis = momentum_basis(sector)
    external = external_wigner_momentum(sector)
    nloops = length(basis) - 1
    iszero(nloops) && return false

    kinematic = kinematic_factor(sector)
    length(kinematic) == 1 || return false
    reference_monomial, reference_coefficient = only(kinematic)
    reference_coefficient == one(ComplexRationals) || return false

    permutation = collect(1:nloops)
    shifts = zeros(Int, nloops)
    signs = ones(Int, nloops)
    for slot in 1:nloops
        signs[slot] = -1
        transform = loop_permutation_transform(basis, external, permutation, signs, shifts)
        transformed_monomial = transform_loop_momenta(monomial, transform)
        if transformed_monomial == monomial
            transformed_support, support_factor = _transform_frequency_support(
                frequency_support(sector), transform
            )
            if transformed_support == frequency_support(sector)
                transformed_kinematic = transform_loop_momenta(kinematic, transform)
                if length(transformed_kinematic) == 1
                    reflected_monomial, kinematic_factor = only(transformed_kinematic)
                    if reflected_monomial == reference_monomial
                        factor =
                            kinematic_factor * convert(ComplexRationals, support_factor)
                        factor != one(ComplexRationals) && return true
                    end
                end
            end
        end
        signs[slot] = 1
    end
    return false
end

function _merge_kernel_polynomial!(
    terms::Dict{CollisionKernelSector{S},OccupationPolynomial{C,S}},
    sector::CollisionKernelSector{S},
    polynomial::OccupationPolynomial{C,S},
) where {C<:Number,S<:Statistics}
    iszero(polynomial) && return terms
    if haskey(terms, sector)
        combined = terms[sector] + polynomial
        if iszero(combined)
            delete!(terms, sector)
        else
            terms[sector] = combined
        end
    else
        terms[sector] = polynomial
    end
    return terms
end

# `quotient_loop_momenta` promotes its output coefficient domain with
# `ComplexRationals`, so this method is the canonical insertion path for quotient
# atoms. Any scalar generated by the chosen loop gauge is moved from the transformed
# kinematic polynomial into the occupation-polynomial coefficient before sector
# identity is compared.
function _push_kernel_polynomial!(
    terms::Dict{CollisionKernelSector{S},OccupationPolynomial{C,S}},
    sector::CollisionKernelSector{S},
    polynomial::OccupationPolynomial{C,S},
) where {C<:Complex,S<:Statistics}
    iszero(polynomial) && return terms

    for (kinematic_monomial, kinematic_coefficient) in kinematic_factor(sector)
        iszero(kinematic_coefficient) && continue
        unit_kinematic = MomentumPolynomial(kinematic_monomial, one(ComplexRationals))
        unit_sector = CollisionKernelSector{S}(
            parameters(sector),
            momentum_basis(sector),
            external_wigner_momentum(sector),
            unit_kinematic,
            frequency_support(sector),
        )
        scale = convert(C, kinematic_coefficient)
        contributions = Pair{OccupationMonomial{S},C}[]
        sizehint!(contributions, length(polynomial))
        for (monomial, coefficient) in polynomial
            _vanishes_under_loop_reflection(unit_sector, monomial) && continue
            push!(contributions, monomial => scale * coefficient)
        end
        isempty(contributions) && continue
        _merge_kernel_polynomial!(
            terms, unit_sector, OccupationPolynomial{C,S}(contributions)
        )
    end

    return terms
end
