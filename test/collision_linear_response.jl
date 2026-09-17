using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields linear_response_ϕ::Boson

struct LinearResponseCustomMoment
    label::Symbol
end

function linear_response_fixture()
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
        KC.OccupationPolynomial(KC.OccupationAtom(linear_response_ϕ, momentum), one(C))
    polynomial = -4 * n(k) * n(q)
    kernel = KC.CollisionKernel{C,Boson,1,0,KC.HomogeneousWignerContext}(
        Dict(sector => polynomial),
        linear_response_ϕ,
        parameter,
        KC.HomogeneousWignerContext(),
    )
    return (; kernel, sector, k, q, n)
end

@testset "exact occupation Fréchet derivative" begin
    fixture = linear_response_fixture()
    nk = fixture.n(fixture.k)
    nq = fixture.n(fixture.q)
    polynomial = 2 * nk * nk * nq + 3 * nq

    linearization = @inferred KC.occupation_linearization(polynomial)
    terms = Dict(KC.occupation_linearization_terms(linearization))
    k_atom = KC.OccupationAtom(linear_response_ϕ, fixture.k)
    q_atom = KC.OccupationAtom(linear_response_ϕ, fixture.q)

    @test length(linearization) == 2
    @test terms[k_atom] == 4 * nk * nq
    @test terms[q_atom] == 2 * nk * nk + 3 * one(typeof(nk))

    constant = 7 * one(typeof(nk))
    @test isempty(@inferred KC.occupation_linearization(constant))
end

@testset "collision kernel linearizes without changing phase-space provenance" begin
    fixture = linear_response_fixture()
    linearized = @inferred KC.linearize_collision_kernel(fixture.kernel)

    @test KC.target_family(linearized) == linear_response_ϕ
    @test KC.parameters(linearized) == KC.parameters(fixture.kernel)
    @test KC.order(linearized) == KC.order(fixture.kernel)
    @test KC.gradient_order(linearized) == KC.gradient_order(fixture.kernel)
    @test KC.wigner_context(linearized) == KC.wigner_context(fixture.kernel)

    response = only(values(KC.linearized_collision_terms(linearized)))
    terms = Dict(KC.occupation_linearization_terms(response))
    k_atom = KC.OccupationAtom(linear_response_ϕ, fixture.k)
    q_atom = KC.OccupationAtom(linear_response_ϕ, fixture.q)

    @test terms[k_atom] == -4 * fixture.n(fixture.q)
    @test terms[q_atom] == -4 * fixture.n(fixture.k)
end

@testset "formal number, energy, and user moment projections" begin
    fixture = linear_response_fixture()

    number_projection = @inferred KC.number_moment_projection(fixture.kernel)
    @test KC.moment_test_function(number_projection) isa KC.NumberMoment
    @test KC.projected_collision(number_projection) === fixture.kernel
    @test KC.moment_weight(number_projection, fixture.sector) == 1 // 1

    energy_projection = @inferred KC.energy_moment_projection(fixture.kernel)
    @test KC.moment_test_function(energy_projection) isa KC.EnergyMoment
    @test KC.moment_weight(energy_projection, fixture.sector) ==
        KC.EnergyForm(KC.DispersionAtom(linear_response_ϕ, fixture.k))

    descriptor = LinearResponseCustomMoment(:quadrupole)
    custom_projection = @inferred KC.project_collision_moment(fixture.kernel, descriptor)
    @test KC.moment_test_function(custom_projection) === descriptor
    @test KC.projected_collision(custom_projection) === fixture.kernel

    linearized_projection = @inferred KC.linearize_collision_moment(energy_projection)
    @test KC.moment_test_function(linearized_projection) isa KC.EnergyMoment
    @test KC.projected_collision(linearized_projection) isa KC.LinearizedCollisionKernel
    @test @inferred(KC.linearize_collision_moment(linearized_projection)) ===
        linearized_projection
end
