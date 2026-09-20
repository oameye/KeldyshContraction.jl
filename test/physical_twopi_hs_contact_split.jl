using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields contact_split_ψ::Boson contact_split_χ::Boson

const ContactSplitCoeff = Complex{Rational{Int}}

function contact_split_hs_interaction()
    ψc = contact_split_ψ[Classical]
    ψq = contact_split_ψ[Quantum]
    χc = contact_split_χ[Classical]
    χq = contact_split_χ[Quantum]
    forward = ψc^2 * bar(χq) + 2 * ψc * ψq * bar(χc) + (1 // 4) * ψq^2 * bar(χq)
    return ChargedInteractionLagrangian(
        -im * (forward + bar(forward)),
        contact_split_ψ => 1,
        contact_split_χ => 2;
        parameter=:h,
    )
end

function contact_split_kind(edge)
    KC.is_keldysh(edge) && return :F
    KC.is_retarded(edge) && return :R
    KC.is_advanced(edge) && return :A
    return error("unexpected propagator type")
end

# The source uses the oppositely oriented \tilde G = G22. Phase A certified that
# a stored KC atomic R line is source \tilde G_A, while stored A is source \tilde G_R.
function contact_split_atomic_kind(edge)
    kind = contact_split_kind(edge)
    kind === :R && return :A
    kind === :A && return :R
    return :F
end

function contact_split_source_pairs(component)
    result = Dict{Tuple{Symbol,Symbol},ContactSplitCoeff}()
    for (diagram, coefficient) in component
        atomic = nothing
        hs = nothing
        for edge in KC.contractions(diagram)
            family = KC.field_family(edge.out)
            if family == contact_split_ψ
                atomic = contact_split_atomic_kind(edge)
            elseif family == contact_split_χ
                hs = contact_split_kind(edge)
            else
                error("unexpected field family")
            end
        end
        atomic === nothing && error("missing atomic line")
        hs === nothing && error("missing HS line")

        # Physical h² = 1/2 and removal of KC's common -i graph phase.
        source_coefficient =
            convert(ContactSplitCoeff, coefficient) * complex(0 // 1, 1 // 2)
        key = (atomic, hs)
        result[key] = get(result, key, zero(ContactSplitCoeff)) + source_coefficient
    end
    return result
end

function contact_split_atomic_factor(kind::Symbol)
    kind === :F && return complex(1 // 1)
    (kind === :R || kind === :A) && return complex(0 // 1, -1 // 1)
    return error("unexpected atomic source component")
end

function contact_split_hs_factors(kind::Symbol, λ, λbar, γ)
    kind === :F && return complex(1 // 1), complex(-γ)
    kind === :R && return complex(0 // 1, -1 // 1), complex(0 // 1, -1 // 1) * λ
    kind === :A && return complex(0 // 1, -1 // 1), complex(0 // 1, -1 // 1) * λbar
    return error("unexpected HS source component")
end

function contact_split_bare_and_regular(pairs, g, γ)
    λ = complex(g, -γ)
    λbar = conj(λ)
    bare_local = Dict{Symbol,ContactSplitCoeff}()
    regular = Dict{Tuple{Symbol,Symbol},ContactSplitCoeff}()

    for ((atomic, hs), coefficient) in pairs
        atomic_factor = contact_split_atomic_factor(atomic)
        regular_factor, bare_factor = contact_split_hs_factors(hs, λ, λbar, γ)

        bare_local[atomic] =
            get(bare_local, atomic, zero(ContactSplitCoeff)) +
            coefficient * atomic_factor * bare_factor
        key = (atomic, hs)
        regular[key] =
            get(regular, key, zero(ContactSplitCoeff)) +
            coefficient * atomic_factor * regular_factor
    end
    return bare_local, regular
end

function contact_split_equal_time(local_terms, F, volume)
    R = complex(0 // 1, -volume / 2)
    A = complex(0 // 1, volume / 2)
    return get(local_terms, :F, zero(ContactSplitCoeff)) * F +
           get(local_terms, :R, zero(ContactSplitCoeff)) * R +
           get(local_terms, :A, zero(ContactSplitCoeff)) * A
end

function contact_split_direct_loss_interaction()
    c = contact_split_ψ[Classical]
    q = contact_split_ψ[Quantum]
    plus = KC.Regularisation.Plus
    minus = KC.Regularisation.Minus
    return im * (
        (1 // 2) * bar(c) * bar(q) * (c(minus) * c(minus) + q(minus) * q(minus)) -
        (1 // 2) * c(plus) * q(plus) * (bar(c) * bar(c) + bar(q) * bar(q)) +
        bar(c) * bar(q) * (c(plus) * q(plus) + c(minus) * q(minus))
    )
end

function contact_split_direct_component_coefficients(component)
    result = Dict{Symbol,ContactSplitCoeff}()
    for (diagram, coefficient) in KC.set_reg_to_zero(component)
        edge = only(KC.contractions(diagram))
        kind = contact_split_kind(edge)
        result[kind] =
            get(result, kind, zero(ContactSplitCoeff)) +
            convert(ContactSplitCoeff, coefficient)
    end
    return result
end

# KC's microscopic first-order quartic-loss oracle is stored in the package's native
# propagator convention. The already-certified source translation used throughout the 2PI
# benchmark is Gᴷ -> -F and Gᴿ/Gᴬ -> -i G_R/G_A. Apply that translation here rather than
# choosing a target qq coefficient directly.
function contact_split_microscopic_qq_source(microscopic, γ)
    scale = convert(ContactSplitCoeff, γ)
    return Dict(
        :F => -scale * microscopic[:F],
        :R => -im * scale * microscopic[:R],
        :A => -im * scale * microscopic[:A],
    )
end

@testset "physical HS bare/regular contact split" begin
    Γ2 = @inferred TwoPIEffectiveAction(contact_split_hs_interaction(), Val(2), Val(3))
    Σ = @inferred SelfEnergy(Γ2, contact_split_ψ)

    advanced_pairs = contact_split_source_pairs(KC.advanced_component(Σ))
    retarded_pairs = contact_split_source_pairs(KC.retarded_component(Σ))
    keldysh_pairs = contact_split_source_pairs(KC.keldysh_component(Σ))

    @test advanced_pairs == Dict((:F, :A) => 2, (:A, :F) => 2)
    @test retarded_pairs == Dict((:F, :R) => 2, (:R, :F) => 2)
    @test keldysh_pairs == Dict((:F, :F) => 2, (:R, :R) => 1 // 2, (:A, :A) => 1 // 2)

    g = 3 // 2
    γ = 2 // 5
    λ = complex(g, -γ)
    λbar = conj(λ)

    advanced_local, advanced_regular = contact_split_bare_and_regular(advanced_pairs, g, γ)
    retarded_local, retarded_regular = contact_split_bare_and_regular(retarded_pairs, g, γ)
    keldysh_local, keldysh_regular = contact_split_bare_and_regular(keldysh_pairs, g, γ)

    @test advanced_local == Dict(:F => -2im * λbar, :A => 2im * γ)
    @test retarded_local == Dict(:F => -2im * λ, :R => 2im * γ)
    @test advanced_regular == Dict((:F, :A) => -2im, (:A, :F) => -2im)
    @test retarded_regular == Dict((:F, :R) => -2im, (:R, :F) => -2im)
    @test keldysh_regular ==
        Dict((:F, :F) => 2, (:R, :R) => -(1 // 2), (:A, :A) => -(1 // 2))

    # Direct D0 substitution is exact for the clq/qcl local channels. At equal time these
    # already reproduce the physical Hartree shift Σ0 and loss rate Γ without any extra graph.
    F = 7 // 5
    volume = 1 // 1
    Σ0 = 2g * F
    Γ = 2γ * (F - volume / 2)
    @test contact_split_equal_time(advanced_local, F, volume) == -im * Σ0 + Γ
    @test contact_split_equal_time(retarded_local, F, volume) == -im * Σ0 - Γ

    # In the qq channel a naive singular substitution into the nonlocal GD formula supplies
    # only half of the causal equal-time contact.
    @test keldysh_local == Dict(:F => -2γ, :R => -(1 // 2) * λ, :A => -(1 // 2) * λbar)
    naive_equal_time = contact_split_equal_time(keldysh_local, F, volume)
    @test naive_equal_time == -2γ * F + γ * volume / 2
    @test naive_equal_time != -Γ

    # Independent microscopic oracle: the established finite-Trotter quartic loss vertex
    # contains the full causal contact pair. Translate its generated Keldysh self-energy to
    # the same source qq convention and normalize it with the physical γ. This determines
    # the missing contact rather than inserting a target coefficient by hand.
    loss = InteractionLagrangian(contact_split_direct_loss_interaction(), :γ)
    G_loss = DressedPropagator(
        loss, Val(1), Val(3); simplify=false, preserve_regularisation=true
    )
    Σ_loss = SelfEnergy(G_loss)
    microscopic = contact_split_direct_component_coefficients(KC.keldysh_component(Σ_loss))
    @test microscopic == Dict(:F => 2, :R => -1, :A => 1)

    microscopic_source = contact_split_microscopic_qq_source(microscopic, γ)
    @test microscopic_source ==
        Dict(:F => complex(-2γ), :R => complex(0 // 1, γ), :A => complex(0 // 1, -γ))
    microscopic_equal_time = contact_split_equal_time(microscopic_source, F, volume)
    @test microscopic_equal_time == -2γ * F + γ * volume
    @test microscopic_equal_time == -Γ
    @test microscopic_equal_time - naive_equal_time == γ * volume / 2

    # The source note's printed λbar*G_R + λ*G_A ordering has the opposite causal sign
    # and therefore fails the same microscopic equal-time oracle.
    printed_note_local = Dict(:F => -2γ, :R => -λbar, :A => -λ)
    printed_equal_time = contact_split_equal_time(printed_note_local, F, volume)
    @test printed_equal_time == -2γ * F - γ * volume
    @test printed_equal_time != microscopic_equal_time

    # At unit loss the comparison exposes the contact completion directly: the statistical
    # coefficient is already correct in the naive HS split, while each causal coefficient
    # carries exactly half of the microscopic finite-Trotter weight.
    pure_naive, _ = contact_split_bare_and_regular(keldysh_pairs, 0 // 1, 1 // 1)
    microscopic_unit = contact_split_microscopic_qq_source(microscopic, 1 // 1)
    @test pure_naive == Dict(
        :F => complex(-2 // 1),
        :R => complex(0 // 1, 1 // 2),
        :A => complex(0 // 1, -1 // 2),
    )
    @test microscopic_unit == Dict(
        :F => complex(-2 // 1),
        :R => complex(0 // 1, 1 // 1),
        :A => complex(0 // 1, -1 // 1),
    )
end

@testset "regular HS remainder keeps the Phase-C0 Fourier graph" begin
    Γ2 = TwoPIEffectiveAction(contact_split_hs_interaction(), Val(2), Val(3))
    ΣF = KC.twopi_fourier_self_energy(Γ2, contact_split_ψ)

    for component in
        (KC.keldysh_component(ΣF), KC.retarded_component(ΣF), KC.advanced_component(ΣF))
        for (graph, contributions) in component
            @test graph.external_count == 1
            @test graph.loop_count == 1
            families = [
                KC.field_family(edge.out) for edge in KC.contractions(graph.coordinate)
            ]
            @test count(==(contact_split_ψ), families) == 1
            @test count(==(contact_split_χ), families) == 1
            @test !isempty(contributions)
        end
    end
end
