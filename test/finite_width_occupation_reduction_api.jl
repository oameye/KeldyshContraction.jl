using KeldyshContraction, Test
import KeldyshContraction as KC

@testset "finite-width occupation API remains qualified-public" begin
    names = (
        :FiniteWidthOccupationCollision,
        :finite_width_occupation_terms,
        :finite_width_unsupported_offset_terms,
        :finite_width_unsupported_distribution_terms,
        :finite_width_occupation_collision,
    )
    for name in names
        if VERSION >= v"1.11"
            @test !Base.isexported(KC, name)
            @test Base.ispublic(KC, name)
        else
            @test !(name in Base.names(KC; all=false, imported=false))
        end
    end
end
