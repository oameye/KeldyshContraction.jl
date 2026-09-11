using KeldyshContraction, Test
import KeldyshContraction as KC
using KeldyshContraction:
    Bulk,
    Contraction,
    Diagram,
    LinearMomentum,
    MomentumComponent,
    MomentumMonomial,
    MomentumPolynomial,
    Out,
    In

function fourier_recursively_concrete(@nospecialize(T::Type), seen=Set{Type}())
    isconcretetype(T) || return false
    T in seen && return true
    push!(seen, T)
    if T <: AbstractArray
        fourier_recursively_concrete(eltype(T), seen) || return false
    elseif T <: AbstractDict
        fourier_recursively_concrete(keytype(T), seen) || return false
        fourier_recursively_concrete(valtype(T), seen) || return false
    end
    for FT in fieldtypes(T)
        fourier_recursively_concrete(FT, seen) || return false
    end
    return true
end

function exact_monomial(momentums::Vector{LinearMomentum}, axis::Symbol)
    return MomentumMonomial(
        MomentumComponent[MomentumComponent(p, axis) for p in momentums]
    )
end

function exact_polynomial(momentums::Vector{LinearMomentum}, axis::Symbol)
    return MomentumPolynomial(exact_monomial(momentums, axis), one(KC.ComplexRationals))
end

@qfields fourier_collection_ϕ::Boson
const fourier_collection_c = fourier_collection_ϕ[Classical]
const fourier_collection_q = fourier_collection_ϕ[Quantum]

function aggregation_diagram(; derivative_endpoint::Symbol)
    out_field = if derivative_endpoint === :out
        partial(fourier_collection_c, :x)
    else
        fourier_collection_c
    end
    in_field = if derivative_endpoint === :in
        partial(bar(fourier_collection_q), :x)
    else
        bar(fourier_collection_q)
    end
    contractions = Contraction{Boson}[
        Contraction(out_field(Out()), bar(fourier_collection_q)(Bulk(1))),
        Contraction(fourier_collection_c(Bulk(1)), bar(fourier_collection_q)(Bulk(1))),
        Contraction(fourier_collection_c(Bulk(1)), in_field(In())),
    ]
    return Diagram(contractions, Val(3), Val(0))
end

function parallel_aggregation_diagram(; derivative_endpoint::Symbol, reverse_parallel=false)
    parallel_out = if derivative_endpoint === :out
        partial(fourier_collection_c, :x)(Bulk(2))
    else
        fourier_collection_c(Bulk(2))
    end
    parallel_in = if derivative_endpoint === :in
        partial(bar(fourier_collection_q), :x)(Bulk(1))
    else
        bar(fourier_collection_q)(Bulk(1))
    end
    decorated = Contraction(parallel_out, parallel_in)
    plain = Contraction(fourier_collection_c(Bulk(2)), bar(fourier_collection_q)(Bulk(1)))
    parallel = reverse_parallel ? (plain, decorated) : (decorated, plain)

    contractions = Contraction{Boson}[
        Contraction(fourier_collection_c(Out()), bar(fourier_collection_q)(Bulk(1))),
        parallel[1],
        parallel[2],
        Contraction(fourier_collection_c(Bulk(1)), bar(fourier_collection_q)(Bulk(2))),
        Contraction(fourier_collection_c(Bulk(2)), bar(fourier_collection_q)(In())),
    ]
    return Diagram(contractions, Val(5), Val(1))
end

@testset "derivative-consumed Fourier collection identity" begin
    out_derivative = aggregation_diagram(; derivative_endpoint=:out)
    in_derivative = aggregation_diagram(; derivative_endpoint=:in)
    @test out_derivative != in_derivative

    diagrams = KC.Diagrams{KC.ComplexRationals,Boson,3,0}()
    push!(diagrams, out_derivative, one(KC.ComplexRationals))
    push!(diagrams, in_derivative, one(KC.ComplexRationals))

    transformed = @inferred fourier_transform(diagrams)
    @test transformed isa FourierDiagrams{KC.ComplexRationals,Boson,3,0}
    @test fourier_recursively_concrete(typeof(transformed))
    @test length(transformed) == 1

    graph, contributions = only(transformed)
    @test length(contributions) == 2
    @test all(
        isempty(derivatives(field)) for edge in KC.contractions(graph.coordinate) for
        field in KC.fields(edge)
    )

    external = LinearMomentum([1, 0])
    expected_out = MomentumPolynomial(
        exact_monomial(LinearMomentum[external], :x), complex(0 // 1, 1 // 1)
    )
    expected_in = MomentumPolynomial(
        exact_monomial(LinearMomentum[external], :x), complex(0 // 1, -1 // 1)
    )
    @test Set(contribution.kinematic for contribution in contributions) ==
        Set((expected_out, expected_in))
end

@testset "derivative-exposed parallel-edge symmetry" begin
    out_derivative = parallel_aggregation_diagram(; derivative_endpoint=:out)
    out_reversed = parallel_aggregation_diagram(;
        derivative_endpoint=:out, reverse_parallel=true
    )
    in_derivative = parallel_aggregation_diagram(; derivative_endpoint=:in)

    @test @inferred(fourier_transform(out_derivative)) ==
        @inferred(fourier_transform(out_reversed))

    diagrams = KC.Diagrams{KC.ComplexRationals,Boson,5,1}()
    push!(diagrams, out_derivative, one(KC.ComplexRationals))
    push!(diagrams, in_derivative, one(KC.ComplexRationals))
    transformed = @inferred fourier_transform(diagrams)

    @test length(transformed) == 1
    graph, contributions = only(transformed)
    @test length(contributions) == 2
    @test KC.external_momentum_count(graph) == 1
    @test KC.loop_momentum_count(graph) == 2
    @test fourier_recursively_concrete(typeof(transformed))
end

@qfields fourier_pwave_ψ::Fermion

function single_fourier_contribution(collection)
    @test length(collection) == 1
    graph, contributions = only(collection)
    @test length(contributions) == 1
    return graph, only(contributions)
end

@testset "p-wave Fourier transform before self-energy amputation" begin
    ψ₁ = fourier_pwave_ψ[One]
    ψ₂ = fourier_pwave_ψ[Two]
    ∂xψ₂ = partial(ψ₂, :x)
    vertex = @inferred ψ₁ * ∂xψ₂ * bar(ψ₁) * bar(∂xψ₂)

    L = @inferred InteractionLagrangian(vertex, :γ)
    G = @inferred DressedPropagator(L, Val(1), Val(3); simplify=false)
    Gk = @inferred fourier_transform(G)
    Σk = @inferred SelfEnergy(Gk)

    @test Gk isa FourierDressedPropagator
    @test Σk isa FourierSelfEnergy
    @test fourier_recursively_concrete(typeof(Gk))
    @test fourier_recursively_concrete(typeof(Σk))

    k = LinearMomentum([1, 0])
    q = LinearMomentum([0, 1])

    retarded_graph, retarded = single_fourier_contribution(Σk.retarded)
    keldysh_graph, keldysh = single_fourier_contribution(Σk.keldysh)
    advanced_graph, advanced = single_fourier_contribution(Σk.advanced)

    @test retarded.kinematic == exact_polynomial(LinearMomentum[q, q], :x)
    @test keldysh.kinematic == exact_polynomial(LinearMomentum[k, q], :x)
    @test advanced.kinematic == exact_polynomial(LinearMomentum[k, k], :x)

    @test retarded.coefficient == complex(0 // 1, -1 // 1)
    @test keldysh.coefficient == complex(0 // 1, 1 // 1)
    @test advanced.coefficient == complex(0 // 1, -1 // 1)

    for graph in (retarded_graph, keldysh_graph, advanced_graph)
        @test KC.external_momentum_count(graph) == 1
        @test KC.loop_momentum_count(graph) == 1
        @test length(KC.momentum_basis(graph)) == 2
        @test length(KC.edge_momenta(graph)) == 1
        @test KC.edge_momenta(graph)[1] == q
    end

    coordinate_Σ = @inferred SelfEnergy(G)
    @test_throws ArgumentError fourier_transform(coordinate_Σ)
end
