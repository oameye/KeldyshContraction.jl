"""One physical ψ kinetic monomial multiplied by one opaque regular pair-response component."""
struct ResponseAwareKineticTerm{E1,E2}
    physical::KineticTerm{Boson,E1,E2}
    response_component::PropagatorType.T
    response_momentum::LinearMomentum
end

kinetic_lines(term::ResponseAwareKineticTerm) = kinetic_lines(term.physical)
kinematic_factor(term::ResponseAwareKineticTerm) = kinematic_factor(term.physical)
momentum_basis(term::ResponseAwareKineticTerm) = momentum_basis(term.physical)
function external_wigner_momentum(term::ResponseAwareKineticTerm)
    return external_wigner_momentum(term.physical)
end
response_component(term::ResponseAwareKineticTerm) = term.response_component
response_momentum(term::ResponseAwareKineticTerm) = term.response_momentum

function Base.isequal(
    a::ResponseAwareKineticTerm{E1,E2}, b::ResponseAwareKineticTerm{E1,E2}
) where {E1,E2}
    return isequal(a.physical, b.physical) &&
           a.response_component === b.response_component &&
           isequal(a.response_momentum, b.response_momentum)
end
Base.:(==)(a::ResponseAwareKineticTerm, b::ResponseAwareKineticTerm) = isequal(a, b)
function Base.hash(term::ResponseAwareKineticTerm, h::UInt)
    return hash(
        ResponseAwareKineticTerm,
        hash(term.response_momentum, hash(term.response_component, hash(term.physical, h))),
    )
end

"""Exact physical-line kinetic expression with an unresolved nested `Dreg[Ω[G]]` factor."""
struct ResponseAwareKineticExpression{C<:Number,E1,E2,Ctx<:AbstractWignerContext}
    terms::Dict{ResponseAwareKineticTerm{E1,E2},C}
    context::Ctx
end

function ResponseAwareKineticExpression{C,E1,E2,Ctx}(
    context::Ctx
) where {C<:Number,E1,E2,Ctx<:AbstractWignerContext}
    K = ResponseAwareKineticTerm{E1,E2}
    return ResponseAwareKineticExpression{C,E1,E2,Ctx}(Dict{K,C}(), context)
end

statistics(::ResponseAwareKineticExpression) = Boson
gradient_order(::ResponseAwareKineticExpression) = Val(0)
wigner_context(expression::ResponseAwareKineticExpression) = expression.context
Base.length(expression::ResponseAwareKineticExpression) = length(expression.terms)
Base.isempty(expression::ResponseAwareKineticExpression) = isempty(expression.terms)
Base.iszero(expression::ResponseAwareKineticExpression) = isempty(expression.terms)
Base.iterate(expression::ResponseAwareKineticExpression) = iterate(expression.terms)
function Base.iterate(expression::ResponseAwareKineticExpression, state)
    return iterate(expression.terms, state)
end

function Base.push!(
    expression::ResponseAwareKineticExpression{C,E1,E2,Ctx},
    term::ResponseAwareKineticTerm{E1,E2},
    coefficient::Number,
) where {C<:Number,E1,E2,Ctx<:AbstractWignerContext}
    value = _simplify(convert(C, coefficient))
    iszero(value) && return expression
    if haskey(expression.terms, term)
        combined = _simplify(expression.terms[term] + value)
        if iszero(combined)
            delete!(expression.terms, term)
        else
            expression.terms[term] = combined
        end
    else
        expression.terms[term] = value
    end
    return expression
end

function Base.:+(
    a::ResponseAwareKineticExpression{C1,E1,E2,Ctx},
    b::ResponseAwareKineticExpression{C2,E1,E2,Ctx},
) where {C1<:Number,C2<:Number,E1,E2,Ctx<:AbstractWignerContext}
    isequal(a.context, b.context) || throw(
        ArgumentError("cannot add response-aware expressions with different contexts")
    )
    D = promote_type(C1, C2)
    out = ResponseAwareKineticExpression{D,E1,E2,Ctx}(a.context)
    for (term, coefficient) in a
        push!(out, term, convert(D, coefficient))
    end
    for (term, coefficient) in b
        push!(out, term, convert(D, coefficient))
    end
    return out
end

function Base.:*(
    prefactor::P, expression::ResponseAwareKineticExpression{C,E1,E2,Ctx}
) where {P<:Number,C<:Number,E1,E2,Ctx<:AbstractWignerContext}
    D = promote_type(P, C)
    out = ResponseAwareKineticExpression{D,E1,E2,Ctx}(expression.context)
    p = convert(D, prefactor)
    for (term, coefficient) in expression
        push!(out, term, p * convert(D, coefficient))
    end
    return out
end
function Base.:*(expression::ResponseAwareKineticExpression, prefactor::Number)
    return prefactor * expression
end
function Base.:-(expression::ResponseAwareKineticExpression{C}) where {C<:Number}
    return -one(C) * expression
end
Base.:-(a::ResponseAwareKineticExpression, b::ResponseAwareKineticExpression) = a + (-b)

function _lower_response_aware_kinetic_term(
    graph::ResponseAwareWignerDiagram{E1,E2,Ctx},
    contribution::WignerContribution{C},
    ::Type{D},
) where {C<:Number,D<:Number,E1,E2,Ctx<:AbstractWignerContext}
    response_index = response_edge_index(graph.composite)
    source_edges = contractions(coordinate_diagram(graph))
    routed_momenta = edge_momenta(graph)
    physical_count = E1 - 1
    lines = Vector{KineticLine{Boson}}(undef, physical_count)
    keldysh_count = Int8(0)
    destination = 1

    @inbounds for source in eachindex(source_edges)
        source == response_index && continue
        line, is_statistical = _lower_kinetic_line(
            source_edges[source], routed_momenta[source]
        )
        lines[destination] = line
        destination += 1
        keldysh_count += is_statistical
    end

    physical = KineticTerm{Boson,physical_count,E2}(
        KineticMonomial(lines, Val(physical_count)),
        topology(coordinate_diagram(graph)),
        momentum_basis(graph),
        external_wigner_momentum(graph),
        contribution.kinematic,
    )
    term = ResponseAwareKineticTerm{physical_count,E2}(
        physical, response_component(graph), response_momentum(graph)
    )
    phase = convert(D, (-im)^Int(keldysh_count))
    coefficient = _simplify(convert(D, contribution.coefficient) * phase)
    return term, coefficient
end

function kinetic_expression(
    diagrams::ResponseAwareWignerDiagrams{C,E1,E2,Ctx}
) where {C<:Number,E1,E2,Ctx<:AbstractWignerContext}
    D = kinetic_coefficient_type(C)
    physical_count = E1 - 1
    out = ResponseAwareKineticExpression{D,physical_count,E2,Ctx}(wigner_context(diagrams))
    for (graph, contributions) in diagrams
        for contribution in contributions
            term, coefficient = _lower_response_aware_kinetic_term(graph, contribution, D)
            push!(out, term, coefficient)
        end
    end
    return out
end

"""
Response-aware physical ψ self-energy after spectral/statistical lowering of physical lines.

`polarization` is the ordinary kinetic lowering of the nested χ-target polarization, whose
internal propagator lines are required to be ψ. It is retained as a nested response input and is
never interpreted as a χ distribution or χ quasiparticle collision.
"""
struct ResponseAwareKineticSelfEnergy{C<:Number,O,E1,E2,R,P,Ctx<:AbstractWignerContext}
    keldysh::ResponseAwareKineticExpression{C,E1,E2,Ctx}
    retarded::ResponseAwareKineticExpression{C,E1,E2,Ctx}
    advanced::ResponseAwareKineticExpression{C,E1,E2,Ctx}
    response::R
    polarization::P
    target::FieldFamily{Boson}
    context::Ctx
end

statistics(::ResponseAwareKineticSelfEnergy) = Boson
order(::ResponseAwareKineticSelfEnergy{C,O}) where {C,O} = O
target_family(self_energy::ResponseAwareKineticSelfEnergy) = self_energy.target
function response_family(self_energy::ResponseAwareKineticSelfEnergy)
    return response_family(self_energy.response)
end
function response_polarization(self_energy::ResponseAwareKineticSelfEnergy)
    return self_energy.polarization
end
function response_coherent_parameter(self_energy::ResponseAwareKineticSelfEnergy)
    return response_coherent_parameter(self_energy.response)
end
function response_loss_parameter(self_energy::ResponseAwareKineticSelfEnergy)
    return response_loss_parameter(self_energy.response)
end
keldysh_component(self_energy::ResponseAwareKineticSelfEnergy) = self_energy.keldysh
retarded_component(self_energy::ResponseAwareKineticSelfEnergy) = self_energy.retarded
advanced_component(self_energy::ResponseAwareKineticSelfEnergy) = self_energy.advanced
gradient_order(::ResponseAwareKineticSelfEnergy) = Val(0)
wigner_context(self_energy::ResponseAwareKineticSelfEnergy) = self_energy.context

function _assert_pure_physical_polarization(
    polarization::KineticSelfEnergy, target::FieldFamily
)
    for expression in (polarization.keldysh, polarization.retarded, polarization.advanced)
        for (term, _) in expression
            for line in kinetic_lines(term)
                isequal(line.family, target) || throw(
                    ArgumentError(
                        "nested HS polarization contains a nonphysical kinetic line"
                    ),
                )
            end
        end
    end
    return polarization
end

function kinetic_expression(
    self_energy::ResponseAwareWignerSelfEnergy{C,O,E1,E2,R,Ctx}
) where {C<:Number,O,E1,E2,R,Ctx<:AbstractWignerContext}
    polarization_wigner = wigner_transform(response_polarization(self_energy.response))
    polarization_kinetic = kinetic_expression(polarization_wigner)
    _assert_pure_physical_polarization(polarization_kinetic, target_family(self_energy))

    keldysh = kinetic_expression(self_energy.keldysh)
    retarded = kinetic_expression(self_energy.retarded)
    advanced = kinetic_expression(self_energy.advanced)
    D = kinetic_coefficient_type(C)
    physical_count = E1 - 1
    return ResponseAwareKineticSelfEnergy{
        D,O,physical_count,E2,R,typeof(polarization_kinetic),Ctx
    }(
        keldysh,
        retarded,
        advanced,
        self_energy.response,
        polarization_kinetic,
        target_family(self_energy),
        self_energy.context,
    )
end

"""Affine homogeneous KB collision with an exact nested pair response."""
struct ResponseAwareOffShellCollisionExpression{
    C<:Number,O,E1,E2,R,P,Ctx<:AbstractWignerContext
}
    offset::ResponseAwareKineticExpression{C,E1,E2,Ctx}
    distribution_coefficient::ResponseAwareKineticExpression{C,E1,E2,Ctx}
    response::R
    polarization::P
    target::FieldFamily{Boson}
    context::Ctx
end

statistics(::ResponseAwareOffShellCollisionExpression) = Boson
order(::ResponseAwareOffShellCollisionExpression{C,O}) where {C,O} = O
target_family(collision::ResponseAwareOffShellCollisionExpression) = collision.target
function response_family(collision::ResponseAwareOffShellCollisionExpression)
    return response_family(collision.response)
end
function response_polarization(collision::ResponseAwareOffShellCollisionExpression)
    return collision.polarization
end
function response_coherent_parameter(collision::ResponseAwareOffShellCollisionExpression)
    return response_coherent_parameter(collision.response)
end
function response_loss_parameter(collision::ResponseAwareOffShellCollisionExpression)
    return response_loss_parameter(collision.response)
end
gradient_order(::ResponseAwareOffShellCollisionExpression) = Val(0)
wigner_context(collision::ResponseAwareOffShellCollisionExpression) = collision.context
collision_offset(collision::ResponseAwareOffShellCollisionExpression) = collision.offset
function collision_distribution_coefficient(
    collision::ResponseAwareOffShellCollisionExpression
)
    return collision.distribution_coefficient
end

function off_shell_collision_expression(
    self_energy::ResponseAwareKineticSelfEnergy{C,O,E1,E2,R,P,Ctx}
) where {C<:Number,O,E1,E2,R,P,Ctx<:AbstractWignerContext}
    D = promote_type(C, ComplexRationals)
    offset = convert(D, im) * self_energy.keldysh
    distribution_coefficient =
        -convert(D, im) * (self_energy.retarded - self_energy.advanced)
    return ResponseAwareOffShellCollisionExpression{D,O,E1,E2,R,P,Ctx}(
        offset,
        distribution_coefficient,
        self_energy.response,
        self_energy.polarization,
        target_family(self_energy),
        self_energy.context,
    )
end
