using KeldyshContraction, Test

@testset "self_energy_type" begin
    using KeldyshContraction: self_energy_type, SmallCollections, In, Out, Keldysh
    dict = SmallCollections.SmallDict{2}(:in => Keldysh, :out => Keldysh)
    @test_throws "should be zero" self_energy_type(dict)
end
