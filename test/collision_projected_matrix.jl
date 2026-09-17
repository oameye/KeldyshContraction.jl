using KeldyshContraction, Test
import KeldyshContraction as KC

@qfields projected_matrix_ϕ::Boson

function projected_matrix_fixture()
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
        KC.OccupationPolynomial(KC.OccupationAtom(projected_matrix_ϕ, momentum), one(C))
    kernel = KC.CollisionKernel{C,Boson,1,0,KC.HomogeneousWignerContext}(
        Dict(sector => -4 * n(k) * n(q)),
        projected_matrix_ϕ,
        parameter,
        KC.HomogeneousWignerContext(),
    )
    k_atom = KC.OccupationAtom(projected_matrix_ϕ, k)
    q_atom = KC.OccupationAtom(projected_matrix_ϕ, q)
    background = KC.OccupationBackground(Rational{Int64}) do atom
        return isequal(atom, k_atom) ? 2 // 1 : 5 // 1
    end
    response = KC.linearize_collision_kernel(kernel, background)
    return (; response, sector, k_atom, q_atom)
end

@testset "explicit left/right projected collision matrix" begin
    fixture = projected_matrix_fixture()
    basis = KC.CollisionProjectionBasis((:number, :energy), (:δk, :δq))

    closure = KC.CollisionPerturbationClosure(Rational{Int64}) do right, atom
        if right === :δk
            return isequal(atom, fixture.k_atom) ? 1 // 1 : 0 // 1
        end
        return isequal(atom, fixture.q_atom) ? 1 // 1 : 0 // 1
    end

    functional = KC.CollisionProjectionFunctional(Rational{Int64}) do left, _, _, response
        return left === :number ? response : 2 * response
    end

    projected = @inferred KC.projected_collision_matrix(
        fixture.response, basis, closure, functional
    )

    @test KC.left_projection_basis(projected) == (:number, :energy)
    @test KC.right_perturbation_basis(projected) == (:δk, :δq)
    @test KC.matrix(projected) == [-20//1 -8//1; -40//1 -16//1]
end

@testset "closure and projection remain explicit" begin
    fixture = projected_matrix_fixture()
    closure = KC.CollisionPerturbationClosure(Rational{Int64}) do right, atom
        return right === :δk && isequal(atom, fixture.k_atom) ? 3 // 1 : 0 // 1
    end
    functional = KC.CollisionProjectionFunctional(
        Rational{Int64}
    ) do left, sector, atom, response
        @test left === :number
        @test sector == fixture.sector
        @test atom == fixture.k_atom
        return response
    end
    basis = KC.CollisionProjectionBasis((:number,), (:δk,))

    projected = @inferred KC.projected_collision_matrix(
        fixture.response, basis, closure, functional
    )
    @test KC.matrix(projected) == reshape(Rational{Int64}[-60 // 1], 1, 1)
end
