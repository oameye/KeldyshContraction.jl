using BenchmarkTools

include(joinpath(@__DIR__, "graphcombinations_loop_equivalence_probe.jl"))
include(joinpath(@__DIR__, "..", "benchmarks", "fermionic_pwave_loss.jl"))

@testset "GC production collision semantic equivalence" begin
    bosonic_collision = benchmark_collision_fixture()
    bosonic_reduced = KC.reduce_frequency_collision(bosonic_collision)
    bosonic_occupation = KC.occupation_reduced_expression(bosonic_reduced)
    bosonic_nauty = KC.quotient_loop_momenta(bosonic_occupation)
    bosonic_gc = gc_quotient_loop_momenta(bosonic_occupation)
    @test gc_requotient_terms(bosonic_nauty) == KC.loop_quotient_terms(bosonic_gc)

    _, _, _, _, _, _, _, _, _, fermionic_occupation, fermionic_nauty =
        benchmark_fermionic_pwave_fixtures()
    fermionic_gc = gc_quotient_loop_momenta(fermionic_occupation)
    @test gc_requotient_terms(fermionic_nauty) == KC.loop_quotient_terms(fermionic_gc)

    println("production collision semantics: bosonic and fermionic GC orbits match Nauty")
end
