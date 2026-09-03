using KeldyshContraction
using Test
using JET

@static if isempty(VERSION.prerelease)
    @testset "JET report_package" begin
        rep = JET.report_package(KeldyshContraction; target_modules=(KeldyshContraction,))
        @show rep
        @test isempty(JET.get_reports(rep))
    end
end
