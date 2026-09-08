using KeldyshContraction, Test
using KeldyshContraction:
    AffineMomentumRouting,
    Bulk,
    Contraction,
    Diagram,
    Edge,
    FixedVector,
    FourierDiagram,
    In,
    LinearMomentum,
    Out,
    amputated_leg_indices,
    basis_momentum,
    contractions,
    coordinate_diagram,
    edge_momenta,
    exact_affine_momentum_routing,
    external_momentum_count,
    index,
    is_bulk,
    loop_momentum_count,
    momentum_basis,
    positions,
    routing_matrix

function fourier_recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        fourier_recursively_concrete(eltype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        fourier_recursively_concrete(FT, seen) || return false
    end
    return true
end

function add_linear!(out::Vector{Rational{Int}}, momentum::LinearMomentum, sign::Int)
    for i in eachindex(out)
        out[i] += sign * momentum[i]
    end
    return out
end

function bulk_conservation(diagram::FourierDiagram, vertex::Int)
    out = zeros(Rational{Int}, length(momentum_basis(diagram)))
    for (edge, momentum) in zip(
        contractions(coordinate_diagram(diagram)), edge_momenta(diagram)
    )
        out_position, in_position = positions(edge)
        if is_bulk(out_position) && index(out_position) == vertex
            add_linear!(out, momentum, 1)
        end
        if is_bulk(in_position) && index(in_position) == vertex
            add_linear!(out, momentum, -1)
        end
    end
    return LinearMomentum(out)
end

@testset "exact affine routing" begin
    incidence = reshape([-1, 1], 2, 1)
    source = reshape([-1, 1], 2, 1)
    routing = @inferred exact_affine_momentum_routing(incidence, source)

    @test routing isa AffineMomentumRouting
    @test external_momentum_count(routing) == 1
    @test loop_momentum_count(routing) == 0
    @test routing_matrix(routing) == reshape(Rational{Int}[1], 1, 1)
    @test incidence * routing_matrix(routing) ==
        Rational{Int}.(source) * reshape(Rational{Int}[1], 1, 1)
    @test fourier_recursively_concrete(typeof(routing))

    inconsistent = reshape([1, 0], 2, 1)
    @test_throws ArgumentError exact_affine_momentum_routing(zeros(Int, 2, 0), inconsistent)
end

@qfields fourier_ϕ::Boson
const fourier_c = fourier_ϕ[Classical]
const fourier_q = fourier_ϕ[Quantum]

function first_order_boson_diagram()
    contractions = Contraction{Boson}[
        Contraction(fourier_c(Out()), bar(fourier_q)(Bulk(1))),
        Contraction(fourier_c(Bulk(1)), bar(fourier_q)(Bulk(1))),
        Contraction(fourier_c(Bulk(1)), bar(fourier_q)(In())),
    ]
    return Diagram(contractions, Val(3), Val(0))
end

function two_loop_boson_diagram()
    contractions = Contraction{Boson}[
        Contraction(fourier_c(Out()), bar(fourier_q)(Bulk(1))),
        Contraction(fourier_c(Bulk(2)), bar(fourier_q)(Bulk(1))),
        Contraction(fourier_c(Bulk(2)), bar(fourier_q)(Bulk(1))),
        Contraction(fourier_c(Bulk(1)), bar(fourier_q)(Bulk(2))),
        Contraction(fourier_c(Bulk(2)), bar(fourier_q)(In())),
    ]
    return Diagram(contractions, Val(5), Val(1))
end

function order_three_boson_diagram()
    contractions = Contraction{Boson}[
        Contraction(fourier_c(Out()), bar(fourier_q)(Bulk(1))),
        Contraction(fourier_c(Bulk(2)), bar(fourier_q)(Bulk(1))),
        Contraction(fourier_c(Bulk(2)), bar(fourier_q)(Bulk(1))),
        Contraction(fourier_c(Bulk(3)), bar(fourier_q)(Bulk(2))),
        Contraction(fourier_c(Bulk(3)), bar(fourier_q)(Bulk(2))),
        Contraction(fourier_c(Bulk(1)), bar(fourier_q)(Bulk(3))),
        Contraction(fourier_c(Bulk(3)), bar(fourier_q)(In())),
    ]
    return Diagram(contractions, Val(7), Val(3))
end

function order_four_boson_diagram()
    contractions = Contraction{Boson}[
        Contraction(fourier_c(Out()), bar(fourier_q)(Bulk(1))),
        Contraction(fourier_c(Bulk(2)), bar(fourier_q)(Bulk(1))),
        Contraction(fourier_c(Bulk(2)), bar(fourier_q)(Bulk(1))),
        Contraction(fourier_c(Bulk(3)), bar(fourier_q)(Bulk(2))),
        Contraction(fourier_c(Bulk(3)), bar(fourier_q)(Bulk(2))),
        Contraction(fourier_c(Bulk(4)), bar(fourier_q)(Bulk(3))),
        Contraction(fourier_c(Bulk(4)), bar(fourier_q)(Bulk(3))),
        Contraction(fourier_c(Bulk(1)), bar(fourier_q)(Bulk(4))),
        Contraction(fourier_c(Bulk(4)), bar(fourier_q)(In())),
    ]
    return Diagram(contractions, Val(9), Val(6))
end

@testset "explicit two-point Fourier diagram" begin
    diagram = first_order_boson_diagram()
    routed = @inferred FourierDiagram(diagram)

    @test routed isa FourierDiagram{Boson,3,0,Nothing}
    @test coordinate_diagram(routed) == diagram
    @test external_momentum_count(routed) == 1
    @test loop_momentum_count(routed) == 1
    @test length(momentum_basis(routed)) == 2
    @test fourier_recursively_concrete(typeof(routed))

    q = basis_momentum(momentum_basis(routed), 1)
    loop = basis_momentum(momentum_basis(routed), 2)
    for (edge, momentum) in zip(contractions(diagram), edge_momenta(routed))
        out_position, in_position = positions(edge)
        if is_bulk(out_position) && is_bulk(in_position)
            @test momentum == loop
        else
            @test momentum == q
        end
    end
    @test iszero(bulk_conservation(routed, 1))

    routed_again = @inferred FourierDiagram(diagram)
    @test routed_again == routed
    @test hash(routed_again) == hash(routed)

    altered = collect(edge_momenta(routed))
    altered[1] = zero(altered[1])
    changed = FourierDiagram{Boson,3,0,Nothing}(
        diagram,
        momentum_basis(routed),
        FixedVector{3,LinearMomentum}(altered),
        Int16(1),
        Int16(1),
        nothing,
    )
    @test changed != routed
end

@testset "multiple loops and higher-order adapters" begin
    two_loop = @inferred FourierDiagram(two_loop_boson_diagram())
    @test external_momentum_count(two_loop) == 1
    @test loop_momentum_count(two_loop) == 2
    @test all(iszero(bulk_conservation(two_loop, vertex)) for vertex in 1:2)

    order_three = @inferred FourierDiagram(order_three_boson_diagram())
    @test loop_momentum_count(order_three) == 3
    @test all(iszero(bulk_conservation(order_three, vertex)) for vertex in 1:3)

    order_four = @inferred FourierDiagram(order_four_boson_diagram())
    @test loop_momentum_count(order_four) == 4
    @test all(iszero(bulk_conservation(order_four, vertex)) for vertex in 1:4)
end

@testset "construction-order determinism" begin
    original = two_loop_boson_diagram()
    reversed_edges = reverse(collect(contractions(original)))
    reconstructed = Diagram(reversed_edges, Val(5), Val(1))

    expected = @inferred FourierDiagram(original)
    actual = @inferred FourierDiagram(reconstructed)
    @test actual == expected
    @test hash(actual) == hash(expected)
end

@testset "amputated self-energy external flow" begin
    interaction = -(
        1 // 2 * (fourier_c^2 + fourier_q^2) * bar(fourier_c) * bar(fourier_q) +
        1 // 2 * fourier_c * fourier_q * (bar(fourier_c)^2 + bar(fourier_q)^2)
    )
    L = InteractionLagrangian(interaction)
    G = DressedPropagator(L, Val(2), Val(5); simplify=false)
    Σ = SelfEnergy(G)
    diagram = first(first(Σ.retarded))
    routed = @inferred FourierDiagram(diagram)

    @test external_momentum_count(routed) == 1
    out_leg, in_leg = amputated_leg_indices(diagram)
    @test out_leg > 0
    @test in_leg > 0

    q = basis_momentum(momentum_basis(routed), 1)
    for vertex in unique((out_leg, in_leg))
        expected = zeros(Rational{Int}, length(momentum_basis(routed)))
        vertex == out_leg && add_linear!(expected, q, 1)
        vertex == in_leg && add_linear!(expected, q, -1)
        @test bulk_conservation(routed, vertex) == LinearMomentum(expected)
    end
end

@qfields fourier_ψ::Fermion
const fourier_fc = fourier_ψ[Classical]
const fourier_fq = fourier_ψ[Quantum]

@testset "statistics-neutral routing" begin
    fermion_contractions = Contraction{Fermion}[
        Contraction(fourier_fc(Out()), bar(fourier_fq)(Bulk(1))),
        Contraction(fourier_fc(Bulk(1)), bar(fourier_fq)(Bulk(1))),
        Contraction(fourier_fc(Bulk(1)), bar(fourier_fq)(In())),
    ]
    fermion = Diagram(fermion_contractions, Val(3), Val(0))

    boson_routing = @inferred FourierDiagram(first_order_boson_diagram())
    fermion_routing = @inferred FourierDiagram(fermion)

    @test momentum_basis(boson_routing) == momentum_basis(fermion_routing)
    @test edge_momenta(boson_routing) == edge_momenta(fermion_routing)
    @test loop_momentum_count(boson_routing) == loop_momentum_count(fermion_routing)
end
