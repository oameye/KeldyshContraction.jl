using KeldyshContraction, Test
import KeldyshContraction as KC

function sd_recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        sd_recursively_concrete(eltype(T), seen) || return false
    elseif T <: AbstractDict
        sd_recursively_concrete(keytype(T), seen) || return false
        sd_recursively_concrete(valtype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        sd_recursively_concrete(FT, seen) || return false
    end
    return true
end

@qfields sd_ϕ::Boson

function sd_self_energy(coefficient; loss=false)
    c = sd_ϕ[Classical]
    q = sd_ϕ[Quantum]
    interaction = if loss
        plus = KC.Regularisation.Plus
        minus = KC.Regularisation.Minus
        coefficient *
        im *
        (
            bar(c) * bar(q) * (c(minus)^2 + q(minus)^2) -
            c(plus) * q(plus) * (bar(c)^2 + bar(q)^2) +
            2 * bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
        )
    else
        -coefficient *
        ((c^2 + q^2) * bar(c) * bar(q) + c * q * (bar(c)^2 + bar(q)^2))
    end
    L = InteractionLagrangian(interaction, loss ? :γ : :g)
    G = DressedPropagator(L, Val(2), Val(5); simplify=false)
    ΣW = wigner_transform(SelfEnergy(fourier_transform(G)); gradient_order=Val(0))
    return kinetic_expression(ΣW)
end

@testset "exact causal-line decomposition convention" begin
    KΣ = sd_self_energy(1 // 2)
    seed = first(kinetic_lines(first(keys(KΣ.keldysh.terms))))
    C = valtype(typeof(KΣ.retarded.terms))
    half_i = convert(C, (1 // 2) * im)

    function variant(kind, statistical=NoStatisticalWeight)
        return KineticLine{Boson}(
            seed.family,
            seed.out_position,
            seed.in_position,
            seed.regularisation_shift,
            KC.momentum(seed),
            kind,
            statistical,
        )
    end

    rparts = @inferred KC._spectral_dispersive_components(variant(KineticRetarded), C)
    aparts = @inferred KC._spectral_dispersive_components(variant(KineticAdvanced), C)
    kparts = @inferred KC._spectral_dispersive_components(
        variant(KineticSpectral, DistributionWeight), C
    )

    @test rparts == [CollisionDispersive => one(C), CollisionSpectral => -half_i]
    @test aparts == [CollisionDispersive => one(C), CollisionSpectral => half_i]
    @test kparts == [CollisionSpectral => one(C)]
end

@testset "i(R-A) cancels dispersive parts exactly" begin
    KΣ = sd_self_energy(1 // 2)
    source = first(keys(KΣ.keldysh.terms))
    source_lines = collect(kinetic_lines(source))
    seed = first(source_lines)

    function causal_variant(kind)
        lines = copy(source_lines)
        lines[1] = KineticLine{Boson}(
            seed.family,
            seed.out_position,
            seed.in_position,
            seed.regularisation_shift,
            KC.momentum(seed),
            kind,
            NoStatisticalWeight,
        )
        return KineticTerm(
            KineticMonomial(lines, Val(3)),
            source.topology,
            source.basis,
            source.external_momentum,
            source.kinematic,
        )
    end

    C = valtype(typeof(KΣ.retarded.terms))
    expression = typeof(KΣ.retarded)(wigner_context(KΣ))
    push!(expression, causal_variant(KineticRetarded), convert(C, im))
    push!(expression, causal_variant(KineticAdvanced), convert(C, -im))

    decomposed = @inferred KC._spectral_dispersive_expression(expression)
    @test length(decomposed) == 1
    term, coefficient = only(decomposed)
    @test coefficient == one(C)
    @test all(==(CollisionSpectral), spectral_dispersive_kinds(term))
end

@testset "complete collision decomposes before shell projection" begin
    KΣ = sd_self_energy(1 // 2)
    collision = @inferred off_shell_collision_expression(KΣ)
    decomposed = @inferred spectral_dispersive_collision(collision)

    @test decomposed isa SpectralDispersiveCollision{
        KC.ComplexRationals,Boson,2,3,1,0,HomogeneousWignerContext
    }
    @test KC.statistics(decomposed) === Boson
    @test KC.order(decomposed) == 2
    @test @inferred(target_family(decomposed)) === sd_ϕ
    @test @inferred(parameters(decomposed)) == parameters(collision)
    @test @inferred(gradient_order(decomposed)) == Val(0)
    @test @inferred(wigner_context(decomposed)) == wigner_context(collision)
    @test sd_recursively_concrete(typeof(decomposed))

    expressions = (collision_offset(decomposed), collision_distribution_coefficient(decomposed))
    @test all(expression -> !iszero(expression), expressions)
    @test all(
        expression -> all(
            term -> all(
                line -> kinetic_line_kind(line) === KineticSpectral,
                kinetic_lines(term.carrier),
            ),
            keys(expression.terms),
        ),
        expressions,
    )
    kinds = Set(
        kind for expression in expressions for (term, _) in expression for
        kind in spectral_dispersive_kinds(term)
    )
    @test CollisionDispersive in kinds
    @test CollisionSpectral in kinds

    again = @inferred spectral_dispersive_collision(collision)
    @test decomposed == again
    @test isequal(decomposed, again)
    @test hash(decomposed) == hash(again)
end

@testset "coefficient domain and regularisation provenance are preserved" begin
    exact = @inferred spectral_dispersive_collision(
        off_shell_collision_expression(sd_self_energy(1 // 2))
    )
    floating = @inferred spectral_dispersive_collision(
        off_shell_collision_expression(sd_self_energy(0.5))
    )
    @test valtype(typeof(collision_offset(exact).terms)) === KC.ComplexRationals
    @test valtype(typeof(collision_distribution_coefficient(floating).terms)) === ComplexF64

    loss_collision = off_shell_collision_expression(sd_self_energy(1 // 2; loss=true))
    loss_decomposed = @inferred spectral_dispersive_collision(loss_collision)
    source_shifts = Set(
        regularisation_shift(line) for expression in (
            collision_offset(loss_collision), collision_distribution_coefficient(loss_collision)
        ) for (term, _) in expression for line in kinetic_lines(term)
    )
    decomposed_shifts = Set(
        regularisation_shift(line) for expression in (
            collision_offset(loss_decomposed),
            collision_distribution_coefficient(loss_decomposed),
        ) for (term, _) in expression for line in kinetic_lines(term.carrier)
    )
    @test source_shifts == decomposed_shifts
    @test any(x -> !iszero(x), source_shifts)
end

@qfields sd_ψ::Fermion sd_χ::Fermion

@testset "fermionic p-wave kinematics survive exact causal decomposition" begin
    ψ₁ = sd_ψ[One]
    χ₂ = sd_χ[Two]
    ∂xχ₂ = partial(χ₂, :x)
    L = InteractionLagrangian((1 // 2) * ψ₁ * ∂xχ₂ * bar(ψ₁) * bar(∂xχ₂), :d)
    G = DressedPropagator(L, Val(1), Val(3); target=sd_ψ, simplify=false)
    ΣW = wigner_transform(SelfEnergy(fourier_transform(G)); gradient_order=Val(0))
    KΣ = @inferred kinetic_expression(ΣW)
    collision = @inferred off_shell_collision_expression(KΣ)
    decomposed = @inferred spectral_dispersive_collision(collision)

    @test KC.statistics(decomposed) === Fermion
    @test @inferred(target_family(decomposed)) === sd_ψ
    @test sd_recursively_concrete(typeof(decomposed))

    source_kinematics = Set(
        KC.kinematic_factor(term) for expression in (
            collision_offset(collision), collision_distribution_coefficient(collision)
        ) for (term, _) in expression
    )
    decomposed_kinematics = Set(
        KC.kinematic_factor(term) for expression in (
            collision_offset(decomposed), collision_distribution_coefficient(decomposed)
        ) for (term, _) in expression
    )
    @test decomposed_kinematics ⊆ source_kinematics
    @test any(
        kinematic -> any(pair -> !isempty(first(pair)), kinematic), decomposed_kinematics
    )
end
