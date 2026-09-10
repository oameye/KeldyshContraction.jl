"""
Exact zeroth-gradient Kadanoff--Baym collision expression before quasiparticle projection.

The collision identity is stored affinely in the external statistical distribution,
`I_coll = I₀ + F_target(k) I₁`, with `I₀ = iΣᴷ` and `I₁ = -A_Σ`.
Here `A_Σ = im * (Σᴿ - Σᴬ)`.

The external distribution is deliberately not represented as an internal `KineticLine`: it
belongs to the external leg and is applied only at the later occupation/on-shell reduction
stage. The physical target family is inherited from `KineticSelfEnergy` rather than repeated
at the collision boundary.
"""
struct OffShellCollisionExpression{
    C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext
}
    offset::KineticExpression{C,S,E1,E2,G,Ctx}
    distribution_coefficient::KineticExpression{C,S,E1,E2,G,Ctx}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::OffShellCollisionExpression{C,S,O}) where {C,S,O} = O
statistics(::OffShellCollisionExpression{C,S}) where {C,S} = S
parameters(collision::OffShellCollisionExpression) = collision.parameter
gradient_order(::OffShellCollisionExpression{C,S,O,E1,E2,G}) where {C,S,O,E1,E2,G} = Val(G)
wigner_context(collision::OffShellCollisionExpression) = collision.context
target_family(collision::OffShellCollisionExpression) = collision.target

"""Return the distribution-independent part `I₀ = im * Σᴷ`."""
collision_offset(collision::OffShellCollisionExpression) = collision.offset

"""Return the coefficient `I₁ = -A_Σ` of the external statistical distribution."""
function collision_distribution_coefficient(collision::OffShellCollisionExpression)
    return collision.distribution_coefficient
end

function Base.isequal(
    a::OffShellCollisionExpression{C,S,O,E1,E2,G,Ctx},
    b::OffShellCollisionExpression{C,S,O,E1,E2,G,Ctx},
) where {C,S,O,E1,E2,G,Ctx}
    return isequal(a.offset, b.offset) &&
           isequal(a.distribution_coefficient, b.distribution_coefficient) &&
           isequal(a.target, b.target) &&
           isequal(a.parameter, b.parameter) &&
           isequal(a.context, b.context)
end
Base.:(==)(a::OffShellCollisionExpression, b::OffShellCollisionExpression) = isequal(a, b)
function Base.hash(collision::OffShellCollisionExpression, h::UInt)
    return hash(
        OffShellCollisionExpression,
        hash(
            (
                collision.offset,
                collision.distribution_coefficient,
                collision.target,
                collision.parameter,
                collision.context,
            ),
            h,
        ),
    )
end

"""
    off_shell_collision_expression(Σ::KineticSelfEnergy)

Form the exact homogeneous zeroth-gradient Kadanoff--Baym collision identity
`I_coll = iΣᴷ - i F_target(k) (Σᴿ - Σᴬ) = iΣᴷ - F_target(k) A_Σ`.

The complete retarded-minus-advanced self-energy is formed before the external statistical
factor is introduced. No quasiparticle shell, energy delta function, occupation substitution,
or principal-value reduction is performed here.

The physical target family is the one already selected by `DressedPropagator` and preserved
through amputation, Fourier transformation, Wigner transformation, and kinetic lowering.
"""
function off_shell_collision_expression(
    Σ::KineticSelfEnergy{C,S,O,E1,E2,0,Ctx}
) where {C<:Number,S<:Statistics,O,E1,E2,Ctx<:AbstractWignerContext}
    D = kinetic_coefficient_type(C)
    offset = convert(D, im) * Σ.keldysh
    distribution_coefficient = -spectral_self_energy(Σ)
    return OffShellCollisionExpression{D,S,O,E1,E2,0,Ctx}(
        offset,
        distribution_coefficient,
        target_family(Σ),
        Σ.parameter,
        Σ.context,
    )
end

function off_shell_collision_expression(
    ::KineticSelfEnergy{C,S,O,E1,E2,G,Ctx}
) where {C<:Number,S<:Statistics,O,E1,E2,G,Ctx<:AbstractWignerContext}
    return throw(
        ArgumentError(
            "off-shell collision construction supports only Wigner gradient order zero; got $G",
        ),
    )
end
