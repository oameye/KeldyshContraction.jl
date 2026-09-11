using KeldyshContraction, Test
import KeldyshContraction as KC

function collision_recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        collision_recursively_concrete(eltype(T), seen) || return false
    elseif T <: AbstractDict
        collision_recursively_concrete(keytype(T), seen) || return false
        collision_recursively_concrete(valtype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        collision_recursively_concrete(FT, seen) || return false
    end
    return true
end

@qfields collision_ϕ::Boson

function bosonic_loss_kinetic_self_energy(coefficient)
    c = collision_ϕ[Classical]
    q = collision_ϕ[Quantum]
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus
    interaction =
        coefficient *
        im *
        (
            bar(c) * bar(q) * (c(minus)^2 + q(minus)^2) -
            c(plus) * q(plus) * (bar(c)^2 + bar(q)^2) +
            2 * bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
        )
    L = InteractionLagrangian(interaction, :γ)
    G = DressedPropagator(L, Val(1), Val(3); simplify=false)
    ΣW = wigner_transform(SelfEnergy(fourier_transform(G)); gradient_order=Val(0))
    return kinetic_expression(ΣW)
end

function bosonic_elastic_kinetic_self_energy(coefficient)
    c = collision_ϕ[Classical]
    q = collision_ϕ[Quantum]
    interaction =
        -coefficient * ((c^2 + q^2) * bar(c) * bar(q) + c * q * (bar(c)^2 + bar(q)^2))
    L = InteractionLagrangian(interaction, :g)
    G = DressedPropagator(L, Val(2), Val(5); simplify=false)
    ΣW = wigner_transform(SelfEnergy(fourier_transform(G)); gradient_order=Val(0))
    return kinetic_expression(ΣW)
end

function edgewise_spectral_term(term::KineticTerm{S,E1,E2}) where {S<:KC.Statistics,E1,E2}
    lines = KineticLine{S}[
        KineticLine{S}(
            line.family,
            line.out_position,
            line.in_position,
            line.regularisation_shift,
            KC.momentum(line),
            KineticSpectral,
            line.statistical,
        ) for line in kinetic_lines(term)
    ]
    return KineticTerm{S,E1,E2}(
        KineticMonomial(lines, Val(E1)),
        term.topology,
        term.basis,
        term.external_momentum,
        term.kinematic,
    )
end

function edgewise_spectral_product_surrogate(
    expression::KineticExpression{C,S,E1,E2,G,Ctx}
) where {C<:Number,S<:KC.Statistics,E1,E2,G,Ctx<:KC.AbstractWignerContext}
    out = KineticExpression{C,S,E1,E2,G,Ctx}(wigner_context(expression))
    for (term, coefficient) in expression
        push!(out, edgewise_spectral_term(term), coefficient)
    end
    return out
end

@testset "exact off-shell Kadanoff-Baym collision identity" begin
    KΣ = bosonic_loss_kinetic_self_energy(1 // 2)
    collision = @inferred off_shell_collision_expression(KΣ)

    @test collision isa OffShellCollisionExpression{
        KC.ComplexRationals,Boson,1,1,0,0,HomogeneousWignerContext
    }
    @test KC.statistics(collision) === Boson
    @test KC.order(collision) == 1
    @test @inferred(target_family(KΣ)) === collision_ϕ
    @test @inferred(target_family(collision)) === collision_ϕ
    @test @inferred(parameters(collision)) == parameters(KΣ)
    @test @inferred(gradient_order(collision)) == Val(0)
    @test @inferred(wigner_context(collision)) == wigner_context(KΣ)
    @test collision_recursively_concrete(typeof(collision))

    C = valtype(typeof(KΣ.keldysh.terms))
    @test !iszero(@inferred(collision_offset(collision)))
    @test !iszero(@inferred(collision_distribution_coefficient(collision)))
    @test @inferred(collision_offset(collision)) == convert(C, im) * KΣ.keldysh
    @test @inferred(collision_distribution_coefficient(collision)) ==
        -spectral_self_energy(KΣ)

    collision_again = @inferred off_shell_collision_expression(KΣ)
    @test collision == collision_again
    @test isequal(collision, collision_again)
    @test hash(collision) == hash(collision_again)
end

@testset "collision coefficient domain follows the kinetic self-energy" begin
    exact = @inferred off_shell_collision_expression(
        bosonic_loss_kinetic_self_energy(1 // 2)
    )
    floating = @inferred off_shell_collision_expression(
        bosonic_loss_kinetic_self_energy(0.5)
    )

    @test valtype(typeof(collision_offset(exact).terms)) === KC.ComplexRationals
    @test valtype(typeof(collision_distribution_coefficient(exact).terms)) ===
        KC.ComplexRationals
    @test valtype(typeof(collision_offset(floating).terms)) === ComplexF64
    @test valtype(typeof(collision_distribution_coefficient(floating).terms)) === ComplexF64
end

@testset "complete R/A subtraction is not an edgewise imaginary-part product" begin
    KΣ = bosonic_elastic_kinetic_self_energy(1 // 2)
    AΣ = @inferred spectral_self_energy(KΣ)
    edgewise = @inferred edgewise_spectral_product_surrogate(KΣ.retarded)

    @test !iszero(AΣ)
    @test AΣ != edgewise
    @test all(
        term ->
            all(line -> kinetic_line_kind(line) === KineticSpectral, kinetic_lines(term)),
        keys(edgewise.terms),
    )
    @test any(
        term -> any(
            line -> kinetic_line_kind(line) in (KineticRetarded, KineticAdvanced),
            kinetic_lines(term),
        ),
        keys(AΣ.terms),
    )
end

@testset "R/A/K identity agrees with greater/lesser decomposition" begin
    Σgreater = (2 // 3) + (1 // 5) * im
    Σlesser = (-1 // 7) + (3 // 11) * im
    ΣK = Σgreater + Σlesser
    discontinuity = Σgreater - Σlesser

    function rak_collision(F)
        return im * ΣK - im * F * discontinuity
    end
    function greater_lesser_collision(F)
        return im * ((1 - F) * Σgreater + (1 + F) * Σlesser)
    end

    @test rak_collision(5 // 4) == greater_lesser_collision(5 // 4)

    n = 2 // 5
    FB = statistical_from_occupation(Boson, n)
    FF = statistical_from_occupation(Fermion, n)
    @test rak_collision(FB) == 2im * ((1 + n) * Σlesser - n * Σgreater)
    @test rak_collision(FF) == 2im * (n * Σgreater + (1 - n) * Σlesser)
end

@qfields collision_ψ::Fermion collision_χ::Fermion

@testset "fermionic target and p-wave provenance survive collision construction" begin
    ψ₁ = collision_ψ[One]
    χ₂ = collision_χ[Two]
    ∂xχ₂ = partial(χ₂, :x)
    interaction = (1 // 2) * ψ₁ * ∂xχ₂ * bar(ψ₁) * bar(∂xχ₂)
    L = InteractionLagrangian(interaction, :d)
    G = DressedPropagator(L, Val(1), Val(3); target=collision_ψ, simplify=false)
    ΣW = wigner_transform(SelfEnergy(fourier_transform(G)); gradient_order=Val(0))
    KΣ = @inferred kinetic_expression(ΣW)
    collision = @inferred off_shell_collision_expression(KΣ)

    @test collision isa OffShellCollisionExpression{
        KC.ComplexRationals,Fermion,1,1,0,0,HomogeneousWignerContext
    }
    @test KC.statistics(collision) === Fermion
    @test @inferred(target_family(KΣ)) === collision_ψ
    @test @inferred(target_family(collision)) === collision_ψ
    @test collision_recursively_concrete(typeof(collision))

    synthetic = @inferred typeof(KΣ)(
        KΣ.retarded, KΣ.retarded, KΣ.advanced, KΣ.parameter, target_family(KΣ), KΣ.context
    )
    synthetic_collision = @inferred off_shell_collision_expression(synthetic)
    source_kinematics = Set(
        KC.kinematic_factor(term) for
        expression in (synthetic.keldysh, synthetic.retarded, synthetic.advanced) for
        (term, _) in expression
    )
    collision_kinematics = Set(
        KC.kinematic_factor(term) for expression in (
            collision_offset(synthetic_collision),
            collision_distribution_coefficient(synthetic_collision),
        ) for (term, _) in expression
    )
    @test !isempty(collision_kinematics)
    @test all(kinematic -> kinematic in source_kinematics, collision_kinematics)
end
