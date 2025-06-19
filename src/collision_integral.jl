
struct CollisionIntegral
    terms::Vector{Float64}
end
struct BosonicDistributionTerm
    momenta::Vector{Momenta}
end
function Base.isequal(bds1::BosonicDistributionTerm, bds2::BosonicDistributionTerm)
    bool = isequal(bds1.momenta, bds2.momenta)
    return bool
end
function Base.hash(bds::BosonicDistributionTerm, h::UInt)
    return hash(BosonicDistributionTerm, hash(bds.momenta, h))
end

Base.isempty(bds::BosonicDistributionTerm) = isempty(bds.momenta)

struct BosonicDistributionAdd
    terms::Dict{BosonicDistributionTerm,ComplexF64}
end
function BosonicDistributionAdd()
    return BosonicDistributionAdd(Dict{BosonicDistributionTerm,ComplexF64}())
end
function Base.push!(
    collection::BosonicDistributionAdd, diagram::BosonicDistributionTerm, prefactor::Number
)
    if haskey(collection.terms, diagram)
        collection.terms[diagram] += prefactor
    else
        collection.terms[diagram] = prefactor
    end
    return collection
end

function imaginary_part(d::Diagram, coeff::ComplexF64=complex(1.0))
    bds = Vector{Momenta}()
    for edge in contractions(d)
        edgetype = propagator_type(edge)
        if is_keldysh(edgetype)
            push!(bds, edge.momenta)
        end
        coeff *= is_advanced(edgetype) ? 0.5 : -0.5
    end
    return BosonicDistributionTerm(bds), coeff
end

function imaginary_part(ds::Diagrams)
    topo = topologies(ds)
    dict = Dict{keytype(topo),BosonicDistributionAdd}()
    for (t_, keys) in topo
        bda = BosonicDistributionAdd()
        for key in keys
            coeff = ds.diagrams[key]
            bds, coeff′ = imaginary_part(key, coeff)
            push!(bda, bds, coeff′)
        end
        dict[t_] = bda
    end
    return dict
end
