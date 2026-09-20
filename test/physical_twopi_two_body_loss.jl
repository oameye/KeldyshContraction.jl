using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields physical_twopi_ψ::Boson physical_twopi_χ::Boson

function physical_twopi_hs_interaction()
    ψc = physical_twopi_ψ[Classical]
    ψq = physical_twopi_ψ[Quantum]
    χc = physical_twopi_χ[Classical]
    χq = physical_twopi_χ[Quantum]

    # Direct Keldysh rotation of
    #   ψ₊² χ̄₊ - ψ₋² χ̄₋,
    # with ψ± = ψc ± ψq/2 and χ± = χc ± χq/2.
    # The mixed ψc ψq term therefore couples to χ̄c, not χ̄q.
    forward = ψc^2 * bar(χq) + 2 * ψc * ψq * bar(χc) + (1 // 4) * ψq^2 * bar(χq)
    interaction = -im * (forward + bar(forward))
    return ChargedInteractionLagrangian(
        interaction, physical_twopi_ψ => 1, physical_twopi_χ => 2; parameter=:h
    )
end

# Quadratic action kernel. With the source convention
#   i D₀⁻¹ = δ²S/δχ̄δχ,
# this is precisely i D₀⁻¹ rather than D₀⁻¹ itself.
function physical_twopi_iD0_inverse_kernel(g::T, γ::T) where {T<:Real}
    λ = complex(g, -γ)
    λbar = conj(λ)
    C = typeof(λ)
    return C[
        zero(C) -inv(λbar)
        -inv(λ) complex(zero(T), -γ / (g^2 + γ^2))
    ]
end

function physical_twopi_bare_hs_propagator(g::T, γ::T) where {T<:Real}
    λ = complex(g, -γ)
    λbar = conj(λ)
    C = typeof(λ)
    return C[
        complex(-γ) -im * λ
        -im * λbar zero(C)
    ]
end

function physical_twopi_is_line(contraction, family, component)
    return KC.field_family(contraction.out) == family &&
           KC.keldysh_index(contraction.out) === first(component) &&
           KC.keldysh_index(contraction.in) === last(component)
end

function physical_twopi_vertex_slot(field)
    family_offset = if KC.field_family(field) == physical_twopi_ψ
        0
    elseif KC.field_family(field) == physical_twopi_χ
        4
    else
        error("unexpected field family")
    end
    orientation_offset = KC.is_barred(field) ? 2 : 0
    keldysh_offset = KC.keldysh_index(field) === Classical ? 0 : 1
    return family_offset + orientation_offset + keldysh_offset + 1
end

function physical_twopi_vertex_coefficients(expression)
    C = Complex{Rational{Int}}
    result = Dict{NTuple{8,Int},C}()
    for term in KC.terms(expression)
        counts = zeros(Int, 8)
        for field in KC.fields(term)
            counts[physical_twopi_vertex_slot(field)] += 1
        end
        key = Tuple(counts)
        result[key] = get(result, key, zero(C)) + KC.coefficient(term)
    end
    return result
end

function physical_twopi_edge_slot(edge)
    family = KC.field_family(edge.out)
    offset = if family == physical_twopi_ψ
        0
    elseif family == physical_twopi_χ
        3
    else
        error("unexpected field family")
    end
    kind = if KC.is_keldysh(edge)
        1
    elseif KC.is_retarded(edge)
        2
    elseif KC.is_advanced(edge)
        3
    else
        error("unexpected propagator type")
    end
    return offset + kind
end

function physical_twopi_component_coefficients(component)
    C = Complex{Rational{Int}}
    result = Dict{NTuple{6,Int},C}()
    for (diagram, coefficient) in component
        counts = zeros(Int, 6)
        for edge in KC.contractions(diagram)
            counts[physical_twopi_edge_slot(edge)] += 1
        end
        key = Tuple(counts)
        result[key] = get(result, key, zero(C)) + coefficient
    end
    return result
end

@testset "physical two-body-loss HS bare kernel" begin
    g = 3 // 2
    γ = 2 // 5
    λ = complex(g, -γ)
    iD0inv = @inferred physical_twopi_iD0_inverse_kernel(g, γ)
    D0 = @inferred physical_twopi_bare_hs_propagator(g, γ)

    @test iD0inv[1, 1] == 0
    @test iD0inv[1, 2] == -inv(conj(λ))
    @test iD0inv[2, 1] == -inv(λ)
    @test iD0inv[2, 2] == complex(0 // 1, -γ / (g^2 + γ^2))
    @test iD0inv[1, 2] == conj(iD0inv[2, 1])
    @test iD0inv[1, 1] * iD0inv[2, 2] - iD0inv[1, 2] * iD0inv[2, 1] == -inv(g^2 + γ^2)

    # The bare solution of the source Dyson equations is
    # D₀ = [-γ  -iλ; -iλ̄  0].  This independently fixes the qq sign above.
    @test D0 == [complex(-γ) -im * λ; -im * conj(λ) complex(0 // 1)]
    @test iD0inv * D0 == [im 0; 0 im]

    elastic = physical_twopi_iD0_inverse_kernel(g, zero(γ))
    @test elastic == [complex(0 // 1) complex(-inv(g)); complex(-inv(g)) complex(0 // 1)]
end

@testset "physical Keldysh rotation of the charged cubic vertex" begin
    L = physical_twopi_hs_interaction()
    coefficients = physical_twopi_vertex_coefficients(L.lagrangian)

    # Keys count
    # (ψc, ψq, ψ̄c, ψ̄q, χc, χq, χ̄c, χ̄q).
    @test coefficients == Dict(
        (2, 0, 0, 0, 0, 0, 0, 1) => -im,
        (1, 1, 0, 0, 0, 0, 1, 0) => -2im,
        (0, 2, 0, 0, 0, 0, 0, 1) => -(1 // 4) * im,
        (0, 0, 2, 0, 0, 1, 0, 0) => -im,
        (0, 0, 1, 1, 1, 0, 0, 0) => -2im,
        (0, 0, 0, 2, 0, 1, 0, 0) => -(1 // 4) * im,
    )
end

@testset "physical G²D skeleton and self-energies" begin
    L = physical_twopi_hs_interaction()
    Γ2 = @inferred TwoPIEffectiveAction(L, Val(2), Val(3))

    @test !iszero(Γ2)
    @test KC.order(Γ2) == 2
    @test KC.parameters(Γ2) == KC.parameter_monomial(:h)^2

    terms = KC.twopi_terms(Γ2)
    @test !isempty(terms)
    @test length(unique(KC.twopi_topology(diagram) for diagram in keys(terms))) == 1
    for diagram in keys(terms)
        families = [
            KC.field_family(contraction.out) for
            contraction in KC.twopi_contractions(diagram)
        ]
        @test count(==(physical_twopi_ψ), families) == 2
        @test count(==(physical_twopi_χ), families) == 1
    end

    Σformal = @inferred KC.twopi_self_energy(Γ2, physical_twopi_ψ)
    Ωformal = @inferred KC.twopi_self_energy(Γ2, physical_twopi_χ)
    Σ = @inferred SelfEnergy(Γ2, physical_twopi_ψ)
    Ω = @inferred SelfEnergy(Γ2, physical_twopi_χ)

    @test !iszero(Σformal)
    @test !iszero(Ωformal)
    @test !iszero(KC.retarded_component(Σ))
    @test !iszero(KC.advanced_component(Σ))
    @test !iszero(KC.keldysh_component(Σ))
    @test !iszero(KC.retarded_component(Ω))
    @test !iszero(KC.advanced_component(Ω))
    @test !iszero(KC.keldysh_component(Ω))

    # The source vertex carries h = 1/sqrt(2), hence h² = 1/2. In KC's stored
    # vacuum-diagram convention the representative Gcc² Dqq graph has coefficient
    # -2 h²; multiplying by the established i/2 source-convention map gives -i/2.
    ccqq = [
        coefficient for (diagram, coefficient) in terms if begin
            contractions = KC.twopi_contractions(diagram)
            count(
                contraction -> physical_twopi_is_line(
                    contraction, physical_twopi_ψ, (Classical, Classical)
                ),
                contractions,
            ) == 2 &&
                count(
                    contraction -> physical_twopi_is_line(
                        contraction, physical_twopi_χ, (Quantum, Quantum)
                    ),
                    contractions,
                ) == 1
        end
    ]
    @test ccqq == [-2]
    h2 = 1 // 2
    source_Γ2_prefactor = complex(0 // 1, 1 // 2) * only(ccqq) * h2
    @test source_Γ2_prefactor == complex(0 // 1, -1 // 2)

    # Functional differentiation of G²D must retain the exact 2:1 line-cut
    # multiplicity. After h²=1/2 and removal of KC's common -i graph phase this
    # is the source normalization Σ ~ 2 G D and Ω ~ G².
    ψ_coefficients = collect(values(KC.twopi_self_energy_terms(Σformal)))
    χ_coefficients = collect(values(KC.twopi_self_energy_terms(Ωformal)))
    @test any(c -> c == -4im, ψ_coefficients)
    @test any(c -> c == -2im, χ_coefficients)
    @test (-4im * h2) / (-im) == 2
    @test (-2im * h2) / (-im) == 1

    # Exact Keldysh-component oracle from the analytical 2PI derivation. Keys are
    # (ψK, ψR, ψA, χK, χR, χA). The source Σ for the ψ field contains the oppositely
    # oriented \tilde G = G₂₂ line. KC stores one oriented ψ -> ψ̄ propagator, so
    # \tilde G_A(x,x') maps to a stored G_R(x',x) edge and \tilde G_R(x,x') maps
    # to G_A(x',x). The χ line and the Ω sector keep their direct R/A labels.
    # KC also stores a common -i graph phase and symbolic h²; setting h²=1/2 then
    # reproduces the source coefficients exactly.
    ΣA = physical_twopi_component_coefficients(KC.advanced_component(Σ))
    ΣR = physical_twopi_component_coefficients(KC.retarded_component(Σ))
    ΣK = physical_twopi_component_coefficients(KC.keldysh_component(Σ))
    ΩA = physical_twopi_component_coefficients(KC.advanced_component(Ω))
    ΩR = physical_twopi_component_coefficients(KC.retarded_component(Ω))
    ΩK = physical_twopi_component_coefficients(KC.keldysh_component(Ω))

    @test ΣA == Dict((1, 0, 0, 0, 0, 1) => -4im, (0, 1, 0, 1, 0, 0) => -4im)
    @test ΣR == Dict((1, 0, 0, 0, 1, 0) => -4im, (0, 0, 1, 1, 0, 0) => -4im)
    @test ΣK == Dict(
        (1, 0, 0, 1, 0, 0) => -4im, (0, 0, 1, 0, 1, 0) => -im, (0, 1, 0, 0, 0, 1) => -im
    )
    @test ΩA == Dict((1, 0, 1, 0, 0, 0) => -4im)
    @test ΩR == Dict((1, 1, 0, 0, 0, 0) => -4im)
    @test ΩK == Dict(
        (2, 0, 0, 0, 0, 0) => -2im,
        (0, 2, 0, 0, 0, 0) => -(1 // 2) * im,
        (0, 0, 2, 0, 0, 0) => -(1 // 2) * im,
    )
end
