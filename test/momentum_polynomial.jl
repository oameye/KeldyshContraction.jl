using KeldyshContraction, Test
import KeldyshContraction as KC
using KeldyshContraction:
    Bulk,
    Contraction,
    Diagram,
    FourierDiagram,
    In,
    LinearMomentum,
    MomentumComponent,
    MomentumMonomial,
    MomentumPolynomial,
    Out,
    basis_momentum,
    kinematic_factor,
    lower_fourier_derivatives,
    momentum_basis

function polynomial_recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        polynomial_recursively_concrete(eltype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        polynomial_recursively_concrete(FT, seen) || return false
    end
    return true
end

function momentum_term(momentum::LinearMomentum, axis::Symbol, coefficient)
    component = MomentumComponent(momentum, axis)
    monomial = MomentumMonomial(MomentumComponent[component])
    return MomentumPolynomial(monomial, coefficient)
end

@testset "exact momentum polynomial algebra" begin
    p = LinearMomentum([1, 0])
    q = LinearMomentum([0, 1])
    ipx = @inferred momentum_term(p, :x, complex(0 // 1, 1 // 1))
    iqx = @inferred momentum_term(q, :x, complex(0 // 1, 1 // 1))

    @test polynomial_recursively_concrete(typeof(ipx))
    @test iszero(@inferred(ipx - ipx))

    relative = @inferred ipx - iqx
    @test length(relative) == 2

    relative_squared = @inferred relative * relative
    @test length(relative_squared) == 3

    px = MomentumComponent(p, :x)
    qx = MomentumComponent(q, :x)
    expected = MomentumPolynomial{KC.ComplexRationals}([
        MomentumMonomial(MomentumComponent[px, px]) => complex(-1 // 1, 0 // 1),
        MomentumMonomial(MomentumComponent[px, qx]) => complex(2 // 1, 0 // 1),
        MomentumMonomial(MomentumComponent[qx, qx]) => complex(-1 // 1, 0 // 1),
    ])
    @test relative_squared == expected
    @test hash(relative_squared) == hash(expected)
end

@qfields polynomial_ϕ::Boson
const polynomial_c = polynomial_ϕ[Classical]
const polynomial_q = polynomial_ϕ[Quantum]

function decorate(field, axes)
    result = field
    for axis in axes
        result = partial(result, axis)
    end
    return result
end

function derivative_external_diagram(; out_axes=(), in_axes=())
    out = decorate(polynomial_c, out_axes)
    incoming_bar = decorate(bar(polynomial_q), in_axes)
    contractions = Contraction{Boson}[
        Contraction(out(Out()), bar(polynomial_q)(Bulk(1))),
        Contraction(polynomial_c(Bulk(1)), bar(polynomial_q)(Bulk(1))),
        Contraction(polynomial_c(Bulk(1)), incoming_bar(In())),
    ]
    return Diagram(contractions, Val(3), Val(0))
end

function expected_single_external_factor(routed, axis, coefficient)
    external = basis_momentum(momentum_basis(routed), 1)
    return momentum_term(external, axis, coefficient)
end

@testset "Fourier derivative endpoint signs" begin
    out_x = @inferred FourierDiagram(derivative_external_diagram(; out_axes=(:x,)))
    out_x_lowered = @inferred lower_fourier_derivatives(out_x)
    @test kinematic_factor(out_x_lowered) ==
        expected_single_external_factor(out_x, :x, complex(0 // 1, 1 // 1))

    in_x = @inferred FourierDiagram(derivative_external_diagram(; in_axes=(:x,)))
    in_x_lowered = @inferred lower_fourier_derivatives(in_x)
    @test kinematic_factor(in_x_lowered) ==
        expected_single_external_factor(in_x, :x, complex(0 // 1, -1 // 1))

    out_t = @inferred FourierDiagram(derivative_external_diagram(; out_axes=(:t,)))
    out_t_lowered = @inferred lower_fourier_derivatives(out_t)
    @test kinematic_factor(out_t_lowered) ==
        expected_single_external_factor(out_t, :t, complex(0 // 1, -1 // 1))

    in_t = @inferred FourierDiagram(derivative_external_diagram(; in_axes=(:t,)))
    in_t_lowered = @inferred lower_fourier_derivatives(in_t)
    @test kinematic_factor(in_t_lowered) ==
        expected_single_external_factor(in_t, :t, complex(0 // 1, 1 // 1))
end

@testset "repeated and mixed derivative lowering" begin
    routed = @inferred FourierDiagram(derivative_external_diagram(; out_axes=(:y, :x)))
    lowered = @inferred lower_fourier_derivatives(routed)
    external = basis_momentum(momentum_basis(routed), 1)
    expected_monomial = MomentumMonomial(
        MomentumComponent[MomentumComponent(external, :x), MomentumComponent(external, :y)]
    )
    expected = MomentumPolynomial(expected_monomial, complex(-1 // 1, 0 // 1))

    @test kinematic_factor(lowered) == expected
    @test polynomial_recursively_concrete(typeof(lowered))
end

@qfields polynomial_ψ::Fermion

@testset "fermionic derivative self energy reaches momentum polynomial" begin
    ψ₁ = polynomial_ψ[One]
    ψ₂ = polynomial_ψ[Two]
    ∂xψ₂ = partial(ψ₂, :x)
    vertex = @inferred ψ₁ * ∂xψ₂ * bar(ψ₁) * bar(∂xψ₂)

    L = @inferred InteractionLagrangian(vertex, :γ)
    G = @inferred DressedPropagator(L, Val(1), Val(3); simplify=false)
    Σ = @inferred SelfEnergy(G)

    found_derivative_factor = false
    for diagrams in (Σ.retarded, Σ.keldysh, Σ.advanced)
        for diagram in keys(diagrams.diagrams)
            routed = @inferred FourierDiagram(diagram)
            lowered = @inferred lower_fourier_derivatives(routed)
            polynomial = kinematic_factor(lowered)
            @test polynomial isa MomentumPolynomial{KC.ComplexRationals}
            @test polynomial_recursively_concrete(typeof(lowered))
            found_derivative_factor |= any(
                !isempty(first(term)) for term in polynomial.terms
            )
        end
    end
    @test found_derivative_factor
end
