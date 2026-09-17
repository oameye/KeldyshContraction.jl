include("direct_one_loop_quotient_probe.jl")

function reflected_one_loop_expression(
    expression::KC.OccupationReducedExpression{C,S,O,G,Ctx},
) where {C<:Number,S<:KC.Statistics,O,G,Ctx<:KC.AbstractWignerContext}
    out = Dict{KC.ReducedCollisionSector{S},KC.OccupationPolynomial{C,S}}()
    for (sector, polynomial) in KC.occupation_reduced_terms(expression)
        reflection = one_loop_reflection(sector)
        transformed_sector, support_factor = KC._transform_kernel_sector(sector, reflection)
        reflected_sector = reduced_sector(transformed_sector)
        reflected_terms = Pair{KC.OccupationMonomial{S},C}[]
        sizehint!(reflected_terms, length(polynomial))
        for (monomial, coefficient) in polynomial
            reflected_monomial = KC.transform_loop_momenta(monomial, reflection)
            push!(
                reflected_terms,
                reflected_monomial => coefficient * convert(C, support_factor),
            )
        end
        out[reflected_sector] = KC.OccupationPolynomial{C,S}(reflected_terms)
    end
    return KC.OccupationReducedExpression{C,S,O,G,Ctx}(
        out,
        KC.target_family(expression),
        KC.parameters(expression),
        KC.wigner_context(expression),
    )
end

@testset "direct one-loop literal reflection invariance" begin
    for (label, expression) in all_cases
        reflected = reflected_one_loop_expression(expression)
        @test direct_one_loop_quotient(reflected) == direct_one_loop_quotient(expression)
    end
end
