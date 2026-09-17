using KeldyshContraction, Test
import KeldyshContraction as KC

function spectral_self_energy_is_exported(name)
    return if VERSION >= v"1.11"
        Base.isexported(KC, name)
    else
        name in names(KC; all=false, imported=false)
    end
end

@testset "spectral self-energy kernel public boundary" begin
    for name in (
        :SpectralSelfEnergySector,
        :SpectralSelfEnergyKernel,
        :momentum_basis,
        :external_wigner_momentum,
        :spectral_self_energy_terms,
        :spectral_self_energy_blocked_terms,
        :spectral_self_energy_trotter_terms,
        :spectral_self_energy_kernel,
    )
        @test !spectral_self_energy_is_exported(name)
        VERSION >= v"1.11" && @test Base.ispublic(KC, name)
    end
end
