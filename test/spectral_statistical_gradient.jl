using KeldyshContraction, Test
import KeldyshContraction as KC

@testset "spectral/statistical lowering rejects nonzero Wigner gradient order" begin
    context = HomogeneousWignerContext()
    diagrams = KC.WignerDiagrams{KC.ComplexRationals,Boson,1,0,1,typeof(context)}(context)

    @test gradient_order(diagrams) == Val(1)
    @test_throws ArgumentError kinetic_expression(diagrams)
end
