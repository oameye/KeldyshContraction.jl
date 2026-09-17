using KeldyshContraction, Test
import KeldyshContraction as KC

function projected_matrix_is_exported(name)
    return if VERSION >= v"1.11"
        Base.isexported(KC, name)
    else
        name in names(KC; all=false, imported=false)
    end
end

@testset "projected collision matrix public boundary" begin
    for name in (
        :CollisionProjectionBasis,
        :CollisionPerturbationClosure,
        :CollisionProjectionFunctional,
        :ProjectedCollisionMatrix,
        :left_projection_basis,
        :right_perturbation_basis,
        :perturbation_amplitude,
        :project_collision_channel,
        :projected_collision_matrix,
    )
        @test !projected_matrix_is_exported(name)
        VERSION >= v"1.11" && @test Base.ispublic(KC, name)
    end
end
