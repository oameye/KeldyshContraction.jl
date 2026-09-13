using KeldyshContraction, Test
using KeldyshContraction: Bulk, In, Out, coefficient, fields, terms
import KeldyshContraction as KC

@testset "canonical numeric coefficients" begin
    @qfields family::Boson
    c, q = family[Classical], family[Quantum]

    # IEEE signed zero is not part of the algebraic identity of a symbolic term.
    # Canonical storage must therefore remove it before strict `isequal`/`hash`
    # semantics are applied.
    with_negative_zero = @inferred KC.QMul(complex(-0.0, 1.0), Field{Boson}[c])
    with_positive_zero = @inferred KC.QMul(complex(0.0, 1.0), Field{Boson}[c])
    @test isequal(with_negative_zero, with_positive_zero)
    @test hash(with_negative_zero) == hash(with_positive_zero)
    @test !signbit(real(KC.coefficient(with_negative_zero)))

    # Scalar signed zero canonicalizes to the closed symbolic zero of the same
    # concrete coefficient type.
    scalar_zero = @inferred KC.QMul(-0.0, Field{Boson}[c])
    @test scalar_zero isa KC.QMul{Float64,Boson}
    @test iszero(scalar_zero)
    @test isempty(scalar_zero.args_nc)
    @test !signbit(KC.coefficient(scalar_zero))

    # Diagram coefficients have their own storage boundary and must satisfy the
    # same invariant; otherwise physically identical propagators compare unequal.
    contractions = KC.Contraction{Boson}[
        KC.Contraction(c(Out()), bar(q)),
        KC.Contraction(c(Bulk()), bar(q)(Bulk())),
        KC.Contraction(c(Bulk()), bar(q)(In())),
    ]
    diagram = @inferred KC.Diagram(contractions, Val(3), Val(0))
    diagrams = KC.Diagrams{ComplexF64,Boson,3,0}()
    push!(diagrams, diagram, complex(-0.0, 1.0))
    stored = only(values(diagrams.diagrams))
    @test isequal(stored, complex(0.0, 1.0))
    @test !signbit(real(stored))

    # Direct Dict construction is used by the analytical regression fixtures and
    # must not provide a bypass around the canonical coefficient invariant.
    from_dict = @inferred KC.Diagrams(Dict(diagram => complex(-0.0, 1.0)))
    stored_from_dict = only(values(from_dict.diagrams))
    @test isequal(stored_from_dict, complex(0.0, 1.0))
    @test !signbit(real(stored_from_dict))
end

@testset "public coefficient conversion and accessors" begin
    @qfields family::Boson
    c, q = family[Classical], family[Quantum]

    rational_term = (1 // 2) * c * bar(c)
    float_term = 0.5 * c * bar(c)
    complex_float_term = complex(0.5, 0.25) * c * bar(c)
    sum = float_term + 0.25 * q * bar(q)

    converted = @inferred convert_coefficients(ComplexF64, rational_term)
    @test @inferred(coefficient(converted)) == 0.5 + 0.0im
    @test @inferred(fields(converted)) == fields(rational_term)

    rationalized = @inferred rationalize_coefficients(float_term)
    @test @inferred(coefficient(rationalized)) == 1 // 2
    @test rationalize_coefficients(rational_term) == rational_term

    rationalized_complex = @inferred rationalize_coefficients(complex_float_term)
    @test @inferred(coefficient(rationalized_complex)) == complex(1 // 2, 1 // 4)

    rationalized_sum = @inferred rationalize_coefficients(sum)
    @test @inferred(coefficient(rationalized_sum)) == [1 // 2, 1 // 4]
    @test length(@inferred(terms(rationalized_sum))) == 2
    @test @inferred(terms(rationalized)) == (rationalized,)

    returned_fields = @inferred fields(float_term)
    pop!(returned_fields)
    @test length(fields(float_term)) == 2

    returned_terms = @inferred terms(sum)
    pop!(returned_terms)
    @test length(terms(sum)) == 2
end
