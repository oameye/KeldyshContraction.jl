using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields background_response_ϕ::Boson

function background_response_fixture()
    C = KC.ComplexRationals
    basis = KC.MomentumBasis(2)
    external = basis[1]
    k = KC.basis_momentum(basis, 1)
    q = KC.basis_momentum(basis, 2)
    parameter = KC.ParameterMonomial(:γ)
    kinematic = KC.MomentumPolynomial(KC.MomentumMonomial(), one(C))
    support = KC.FrequencySupport(
        KC.EnergyShell{Boson}[], KC.PrincipalValueSupport{Boson}[]
    )
    sector = KC.CollisionKernelSector{Boson}(parameter, basis, external, kinematic, support)
    n(momentum) =
        KC.OccupationPolynomial(KC.OccupationAtom(background_response_ϕ, momentum), one(C))
    kernel = KC.CollisionKernel{C,Boson,1,0,KC.HomogeneousWignerContext}(
        Dict(sector => -4 * n(k) * n(q)),
        background_response_ϕ,
        parameter,
        KC.HomogeneousWignerContext(),
    )
    return (; kernel, sector, k, q, n)
end

function background_response_state(fixture)
    k_atom = KC.OccupationAtom(background_response_ϕ, fixture.k)
    return KC.OccupationBackground(Rational{Int64}) do atom
        return isequal(atom, k_atom) ? 2 // 1 : 5 // 1
    end
end

@testset "typed occupation background evaluates exact polynomials" begin
    fixture = background_response_fixture()
    background = background_response_state(fixture)
    nk = fixture.n(fixture.k)
    nq = fixture.n(fixture.q)
    k_atom = KC.OccupationAtom(background_response_ϕ, fixture.k)
    q_atom = KC.OccupationAtom(background_response_ϕ, fixture.q)

    @test @inferred(KC.occupation_background_value(background, k_atom)) == 2 // 1
    @test @inferred(KC.occupation_background_value(background, q_atom)) == 5 // 1

    polynomial = 2 * nk * nk * nq + 3 * nq
    @test @inferred(KC.evaluate_occupation_polynomial(polynomial, background)) == 55 // 1

    symbolic = @inferred KC.occupation_linearization(polynomial)
    evaluated = @inferred KC.evaluate_occupation_linearization(symbolic, background)
    terms = Dict(KC.background_occupation_linearization_terms(evaluated))
    @test terms[k_atom] == 40 // 1
    @test terms[q_atom] == 11 // 1

    zero_background = KC.OccupationBackground(Rational{Int64}, _ -> 0 // 1)
    zero_response = @inferred KC.evaluate_occupation_linearization(symbolic, zero_background)
    @test isempty(zero_response)
end

@testset "collision response evaluates around a supplied background" begin
    fixture = background_response_fixture()
    background = background_response_state(fixture)
    k_atom = KC.OccupationAtom(background_response_ϕ, fixture.k)
    q_atom = KC.OccupationAtom(background_response_ϕ, fixture.q)

    symbolic = @inferred KC.linearize_collision_kernel(fixture.kernel)
    evaluated = @inferred KC.evaluate_collision_background(symbolic, background)
    direct = @inferred KC.linearize_collision_kernel(fixture.kernel, background)

    @test KC.target_family(evaluated) == KC.target_family(fixture.kernel)
    @test KC.parameters(evaluated) == KC.parameters(fixture.kernel)
    @test KC.order(evaluated) == KC.order(fixture.kernel)
    @test KC.gradient_order(evaluated) == KC.gradient_order(fixture.kernel)
    @test KC.wigner_context(evaluated) == KC.wigner_context(fixture.kernel)
    @test KC.background_linearized_collision_terms(direct) ==
        KC.background_linearized_collision_terms(evaluated)

    response = only(values(KC.background_linearized_collision_terms(evaluated)))
    terms = Dict(KC.background_occupation_linearization_terms(response))
    @test terms[k_atom] == -20 // 1
    @test terms[q_atom] == -8 // 1

    zero_background = KC.OccupationBackground(Rational{Int64}, _ -> 0 // 1)
    zero_kernel = @inferred KC.linearize_collision_kernel(fixture.kernel, zero_background)
    @test isempty(zero_kernel)
end

@testset "moment projection preserves supplied-background response" begin
    fixture = background_response_fixture()
    background = background_response_state(fixture)

    projection = @inferred KC.energy_moment_projection(fixture.kernel)
    evaluated_projection = @inferred KC.linearize_collision_moment(projection, background)
    @test KC.moment_test_function(evaluated_projection) isa KC.EnergyMoment
    @test KC.projected_collision(evaluated_projection) isa KC.BackgroundLinearizedCollisionKernel
    @test KC.moment_weight(evaluated_projection, fixture.sector) ==
        KC.EnergyForm(KC.DispersionAtom(background_response_ϕ, fixture.k))

    evaluated_kernel = KC.projected_collision(evaluated_projection)
    number_projection = @inferred KC.number_moment_projection(evaluated_kernel)
    @test KC.moment_weight(number_projection, fixture.sector) == 1 // 1
end
