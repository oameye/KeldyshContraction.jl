using KeldyshContraction, Test
import KeldyshContraction as KC

function kinetic_recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        kinetic_recursively_concrete(eltype(T), seen) || return false
    elseif T <: AbstractDict
        kinetic_recursively_concrete(keytype(T), seen) || return false
        kinetic_recursively_concrete(valtype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        kinetic_recursively_concrete(FT, seen) || return false
    end
    return true
end

function assert_single_line_lowering(edge, expected_kind, expected_weight, expected_count)
    routed_momentum = KC.LinearMomentum([1, -1])
    line, count = @inferred KC._lower_kinetic_line(edge, routed_momentum)
    @test kinetic_line_kind(line) === expected_kind
    @test statistical_weight(line) === expected_weight
    @test count == expected_count
    @test KC.momentum(line) == routed_momentum
    return line
end

function kinetic_signature(term)
    retarded = 0
    advanced = 0
    distributed_spectral = 0
    bare_spectral = 0
    for line in kinetic_lines(term)
        kind = kinetic_line_kind(line)
        weight = statistical_weight(line)
        if kind === KineticRetarded
            retarded += 1
        elseif kind === KineticAdvanced
            advanced += 1
        elseif weight === DistributionWeight
            distributed_spectral += 1
        else
            bare_spectral += 1
        end
    end
    return (retarded, advanced, distributed_spectral, bare_spectral)
end

function aggregate_kinetic_signatures(expression, topology_value::Int)
    C = valtype(typeof(expression.terms))
    result = Dict{NTuple{4,Int},C}()
    for (term, coefficient) in expression
        collect(term.topology) == [topology_value] || continue
        signature = kinetic_signature(term)
        result[signature] = get(result, signature, zero(C)) + coefficient
    end
    filter!((pair) -> !iszero(last(pair)), result)
    return result
end

@testset "statistics-dependent occupation convention" begin
    @test statistical_occupation_coefficients(Boson) == (Int8(1), Int8(2))
    @test statistical_occupation_coefficients(Fermion) == (Int8(1), Int8(-2))
    @test @inferred(statistical_from_occupation(Boson, 3 // 5)) == 11 // 5
    @test @inferred(statistical_from_occupation(Fermion, 3 // 5)) == -1 // 5
    @test @inferred(statistical_from_occupation(Boson, 0.25)) == 1.5
    @test @inferred(statistical_from_occupation(Fermion, 0.25)) == 0.5
end

@testset "complete-product imaginary part cannot be reduced edge by edge" begin
    z₁ = 1 + 2im
    z₂ = 3 + 4im
    @test imag(z₁ * z₂) == 10
    @test imag(z₁) * imag(z₂) == 8
    @test imag(z₁ * z₂) != imag(z₁) * imag(z₂)
end

@qfields kinetic_boson_ϕ::Boson

@testset "single-line bosonic R/A/K conversion" begin
    c = kinetic_boson_ϕ[Classical]
    q = kinetic_boson_ϕ[Quantum]
    keldysh = KC.Edge(c(KC.Bulk(1)), bar(c)(KC.Bulk(2)))
    retarded = KC.Edge(c(KC.Bulk(1)), bar(q)(KC.Bulk(2)))
    advanced = KC.Edge(q(KC.Bulk(1)), bar(c)(KC.Bulk(2)))

    kline = assert_single_line_lowering(
        keldysh, KineticSpectral, DistributionWeight, Int8(1)
    )
    rline = assert_single_line_lowering(
        retarded, KineticRetarded, NoStatisticalWeight, Int8(0)
    )
    aline = assert_single_line_lowering(
        advanced, KineticAdvanced, NoStatisticalWeight, Int8(0)
    )

    kr = @inferred kline * rline
    rk = @inferred rline * kline
    @test kr isa KineticMonomial{Boson,2}
    @test kr == rk
    @test @inferred(kr * aline) == @inferred(aline * rk)
end

@testset "equal-time regulator spelling has one kinetic identity" begin
    c = kinetic_boson_ϕ[Classical]
    q = kinetic_boson_ϕ[Quantum]
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus
    y = KC.Bulk(1)
    momentum = KC.LinearMomentum([1])

    retarded_plus_zero = KC.Edge(c(plus)(y), bar(q)(y))
    retarded_zero_minus = KC.Edge(c(y), bar(q)(minus)(y))
    line₁, count₁ = @inferred KC._lower_kinetic_line(retarded_plus_zero, momentum)
    line₂, count₂ = @inferred KC._lower_kinetic_line(retarded_zero_minus, momentum)

    @test count₁ == Int8(0)
    @test count₂ == Int8(0)
    @test regularisation_shift(line₁) == Int8(1)
    @test regularisation_shift(line₂) == Int8(1)
    @test line₁ == line₂
    @test isequal(line₁, line₂)
    @test hash(line₁) == hash(line₂)
end

function elastic_boson_interaction(half)
    c = kinetic_boson_ϕ[Classical]
    q = kinetic_boson_ϕ[Quantum]
    return -(half * (c^2 + q^2) * bar(c) * bar(q) + half * c * q * (bar(c)^2 + bar(q)^2))
end

@testset "bosonic Wigner self-energy lowers with exact Gᴷ = -i F A convention" begin
    interaction = elastic_boson_interaction(1 // 2)
    L = InteractionLagrangian(interaction, :g)
    G = DressedPropagator(L, Val(2), Val(5); simplify=false)
    ΣW = wigner_transform(SelfEnergy(fourier_transform(G)); gradient_order=Val(0))
    KΣ = @inferred kinetic_expression(ΣW)

    @test KΣ isa KineticSelfEnergy
    @test KC.statistics(KΣ) === Boson
    @test KC.order(KΣ) == 2
    @test parameters(KΣ) == parameters(ΣW)
    @test gradient_order(KΣ) == Val(0)
    @test wigner_context(KΣ) == wigner_context(ΣW)
    @test kinetic_recursively_concrete(typeof(KΣ))

    graph, contributions = first(ΣW.keldysh)
    contribution = first(contributions)
    D = KC.kinetic_coefficient_type(typeof(contribution.coefficient))
    term, coefficient = @inferred KC._lower_kinetic_term(graph, contribution, D)
    source_edges = KC.contractions(KC.coordinate_diagram(graph))
    keldysh_count = count(KC.is_keldysh, source_edges)

    @test coefficient ==
        convert(D, contribution.coefficient) * convert(D, (-im)^keldysh_count)
    @test count(
        line -> statistical_weight(line) === DistributionWeight, kinetic_lines(term)
    ) == keldysh_count
    @test all(
        line ->
            statistical_weight(line) !== DistributionWeight ||
            kinetic_line_kind(line) === KineticSpectral,
        kinetic_lines(term),
    )
    @test KC.kinematic_factor(term) == contribution.kinematic
    @test KC.momentum_basis(term) == KC.momentum_basis(graph)
    @test external_wigner_momentum(term) == external_wigner_momentum(graph)

    source_momenta = KC.edge_momenta(graph)
    lowered_momenta = sort!(
        collect(KC.momentum(line) for line in kinetic_lines(term));
        lt=KC._linear_momentum_isless,
    )
    expected_momenta = sort!(collect(source_momenta); lt=KC._linear_momentum_isless)
    @test lowered_momenta == expected_momenta

    discontinuity = @inferred retarded_minus_advanced(KΣ)
    spectral = @inferred spectral_self_energy(KΣ)
    @test spectral == im * discontinuity
    @test iszero(@inferred(KΣ.retarded - KΣ.retarded))

    C = valtype(typeof(KΣ.keldysh.terms))
    unit_i = convert(C, im)
    half_i = convert(C, (1 // 2) * im)
    expected_keldysh = Dict{NTuple{4,Int},C}(
        (0, 0, 3, 0) => -half_i,
        (0, 2, 1, 0) => half_i,
        (2, 0, 1, 0) => half_i,
        (1, 1, 1, 0) => 2 * unit_i,
    )
    @test aggregate_kinetic_signatures(KΣ.keldysh, 3) == expected_keldysh
end

@testset "kinetic coefficient domain follows the Wigner coefficient type" begin
    exact_L = InteractionLagrangian(elastic_boson_interaction(1 // 2), :g)
    exact_G = DressedPropagator(exact_L, Val(1), Val(3); simplify=false)
    exact_ΣW = wigner_transform(
        SelfEnergy(fourier_transform(exact_G)); gradient_order=Val(0)
    )
    exact_KΣ = @inferred kinetic_expression(exact_ΣW)

    float_L = InteractionLagrangian(elastic_boson_interaction(0.5), :g)
    float_G = DressedPropagator(float_L, Val(1), Val(3); simplify=false)
    float_ΣW = wigner_transform(
        SelfEnergy(fourier_transform(float_G)); gradient_order=Val(0)
    )
    float_KΣ = @inferred kinetic_expression(float_ΣW)

    @test valtype(typeof(exact_KΣ.retarded.terms)) === KC.ComplexRationals
    @test valtype(typeof(exact_KΣ.keldysh.terms)) === KC.ComplexRationals
    @test valtype(typeof(float_KΣ.retarded.terms)) === ComplexF64
    @test valtype(typeof(float_KΣ.keldysh.terms)) === ComplexF64
end

@qfields kinetic_loss_ϕ::Boson

@testset "regularised loss tadpole provenance survives the default public path" begin
    c = kinetic_loss_ϕ[Classical]
    q = kinetic_loss_ϕ[Quantum]
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus
    loss =
        im * (
            1 // 2 * bar(c) * bar(q) * (c(minus)^2 + q(minus)^2) -
            1 // 2 * c(plus) * q(plus) * (bar(c)^2 + bar(q)^2) +
            bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
        )
    L = InteractionLagrangian(loss, :γ)
    G = DressedPropagator(L, Val(2), Val(5); simplify=false)
    ΣW = wigner_transform(SelfEnergy(fourier_transform(G)); gradient_order=Val(0))
    KΣ = @inferred kinetic_expression(ΣW)

    @test KC.statistics(KΣ) === Boson
    @test kinetic_recursively_concrete(typeof(KΣ))

    C = valtype(typeof(KΣ.keldysh.terms))
    unit_i = convert(C, im)
    half_i = convert(C, (1 // 2) * im)
    expected_three_line_keldysh = Dict{NTuple{4,Int},C}(
        (0, 0, 3, 0) => -half_i,
        (0, 2, 1, 0) => half_i,
        (2, 0, 1, 0) => half_i,
        (1, 0, 2, 0) => -2 * one(C),
        (0, 1, 2, 0) => 2 * one(C),
        (1, 1, 1, 0) => -2 * unit_i,
    )
    @test aggregate_kinetic_signatures(KΣ.keldysh, 3) == expected_three_line_keldysh

    saw_tadpole = false
    saw_shifted_tadpole_line = false
    for collection in (ΣW.keldysh, ΣW.retarded, ΣW.advanced)
        for (graph, contributions) in collection
            coordinate = KC.coordinate_diagram(graph)
            source_edges = KC.contractions(coordinate)
            is_tadpole = collect(KC.topology(coordinate)) == [2]
            saw_tadpole |= is_tadpole
            D = KC.kinetic_coefficient_type(typeof(first(contributions).coefficient))
            for contribution in contributions
                term, _ = @inferred KC._lower_kinetic_term(graph, contribution, D)
                source_shifts = sort([
                    KC._kinetic_regularisation_shift(edge) for edge in source_edges
                ])
                lowered_shifts = sort([
                    regularisation_shift(line) for line in kinetic_lines(term)
                ])
                @test lowered_shifts == source_shifts

                for line in kinetic_lines(term)
                    iszero(regularisation_shift(line)) && continue
                    @test line.out_position == line.in_position
                    @test kinetic_line_kind(line) in (KineticRetarded, KineticAdvanced)
                    @test statistical_weight(line) === NoStatisticalWeight
                    saw_shifted_tadpole_line |= is_tadpole
                end
            end
        end
    end
    @test saw_tadpole
    @test saw_shifted_tadpole_line
end

@qfields kinetic_fermion_ψ::Fermion

@testset "single-line fermionic R/A/K conversion" begin
    one = kinetic_fermion_ψ[One]
    two = kinetic_fermion_ψ[Two]
    retarded = KC.Edge(one(KC.Bulk(1)), bar(one)(KC.Bulk(2)))
    keldysh = KC.Edge(one(KC.Bulk(1)), bar(two)(KC.Bulk(2)))
    advanced = KC.Edge(two(KC.Bulk(1)), bar(two)(KC.Bulk(2)))

    assert_single_line_lowering(retarded, KineticRetarded, NoStatisticalWeight, Int8(0))
    assert_single_line_lowering(keldysh, KineticSpectral, DistributionWeight, Int8(1))
    assert_single_line_lowering(advanced, KineticAdvanced, NoStatisticalWeight, Int8(0))
end

@testset "fermionic p-wave momentum polynomial survives kinetic lowering" begin
    ψ₁ = kinetic_fermion_ψ[One]
    ψ₂ = kinetic_fermion_ψ[Two]
    ∂xψ₂ = partial(ψ₂, :x)
    interaction = ψ₁ * ∂xψ₂ * bar(ψ₁) * bar(∂xψ₂)
    L = InteractionLagrangian(interaction, :γ)
    G = DressedPropagator(L, Val(1), Val(3); target=kinetic_fermion_ψ, simplify=false)
    ΣW = wigner_transform(SelfEnergy(fourier_transform(G)); gradient_order=Val(0))
    KΣ = @inferred kinetic_expression(ΣW)

    @test KC.statistics(KΣ) === Fermion
    @test KC.order(KΣ) == 1
    @test kinetic_recursively_concrete(typeof(KΣ))

    source_kinematics = Set(
        contribution.kinematic for collection in (ΣW.keldysh, ΣW.retarded, ΣW.advanced) for
        (_, contributions) in collection for contribution in contributions
    )
    lowered_kinematics = Set(
        KC.kinematic_factor(term) for expression in (KΣ.keldysh, KΣ.retarded, KΣ.advanced)
        for (term, _) in expression
    )
    @test lowered_kinematics == source_kinematics
    @test any(
        term -> any(pair -> !isempty(first(pair)), KC.kinematic_factor(term)),
        keys(KΣ.retarded.terms),
    )
end
