using KeldyshContraction, Test

@testset "self_energy_type" begin
    using KeldyshContraction: self_energy_type, OrderedCollections, In, Out, Keldysh
    dict = OrderedCollections.LittleDict(:in => Keldysh, :out => Keldysh)
    @test_throws "should be zero" self_energy_type(dict)
end
