using KeldyshContraction
using Test

import KeldyshContraction as KC

@qfields topology_routing_ϕ::Boson
const topology_routing_c = topology_routing_ϕ[Classical]
const topology_routing_q = topology_routing_ϕ[Quantum]

function elastic_interaction()
    c = topology_routing_c
    q = topology_routing_q
    return -(
        1 // 2 * (c^2 + q^2) * bar(c) * bar(q) + 1 // 2 * c * q * (bar(c)^2 + bar(q)^2)
    )
end

function legacy_external_flow_matrix(diagram)
    incidence = KC.construct_linear_system(diagram.contractions)
    size(incidence, 1) == 2 || error("legacy routing probe expects two bulk vertices")
    outgoing = Int[-1, 0]
    incoming = isodd(first(KC.topology(diagram))) ? Int[0, 1] : Int[1, 0]
    return hcat(outgoing, incidence, incoming)
end

function structural_external_flow_matrix(diagram)
    incidence = KC.construct_linear_system(diagram.contractions)
    out_leg, in_leg = KC.amputated_leg_indices(diagram)
    nvertices = size(incidence, 1)
    1 <= out_leg <= nvertices || error("could not infer outgoing amputated attachment")
    1 <= in_leg <= nvertices || error("could not infer incoming amputated attachment")

    outgoing = zeros(Int, nvertices)
    incoming = zeros(Int, nvertices)
    outgoing[out_leg] = -1
    incoming[in_leg] = 1
    return hcat(outgoing, incidence, incoming), out_leg, in_leg
end

function component_records(name, diagrams)
    records = NamedTuple[]
    for (diagram, _) in diagrams
        legacy = legacy_external_flow_matrix(diagram)
        structural, out_leg, in_leg = structural_external_flow_matrix(diagram)
        topology = collect(KC.topology(diagram))
        push!(
            records,
            (
                component=name,
                topology=topology,
                old_in_leg=isodd(first(topology)) ? 2 : 1,
                out_leg=out_leg,
                in_leg=in_leg,
                matrices_equal=legacy == structural,
            ),
        )
    end
    return records
end

@testset "legacy topology routing is structural external flow" begin
    L = InteractionLagrangian(elastic_interaction())
    G = DressedPropagator(L, Val(2), Val(5); simplify=false)
    Σ = SelfEnergy(G)

    records = NamedTuple[]
    append!(records, component_records(:keldysh, Σ.keldysh))
    append!(records, component_records(:retarded, Σ.retarded))
    append!(records, component_records(:advanced, Σ.advanced))
    isempty(records) && error("self-energy fixture produced no diagrams")

    println("legacy topology routing probe: ", length(records), " diagrams")
    for record in records
        println(
            record.component,
            " topology=",
            record.topology,
            " old attachments=(1,",
            record.old_in_leg,
            ") inferred=(",
            record.out_leg,
            ",",
            record.in_leg,
            ") matrix_equal=",
            record.matrices_equal,
        )
        @test record.out_leg == 1
        @test record.in_leg == record.old_in_leg
        @test record.matrices_equal
    end
end
