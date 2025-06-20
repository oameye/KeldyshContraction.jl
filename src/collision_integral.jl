
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

"""
⃗a ⃗a = - A A -  ⃗r ⃗r
⃗a ⃖r =  A A -  ⃗r ⃖a
"""
function reduce_to_spectral(ds::Diagrams{E1,E2}) where {E1,E2}
    ds′ = Diagrams{E1,E2}()
    for (d, coeff) in ds.diagrams
        contractions = contractions(d)
        types = propagator_type.(contractions)
        directions = direction.(contractions)

        idxs = [1, 2]
        aa_simplification = false
        ar_simplification = false
        while last(idxs) <= E1
            if all(is_advanced, types[idxs]) && all(directions[idxs])
                aa_simplification = true
                break
            elseif types[idxs] == [PropagatorType.Advanced, PropagatorType.Retarded] &&
                directions[idxs] == [true, false]
                ar_simplification = true
                break
            end
            idxs .+= 1
        end
        if last(idxs) > E1
            push!(ds′, d, coeff)
            continue
        elseif aa_simplification
            contractions′ = FixedVector(
                i ∈ idxs ? make_spectral(e) : e for (i, e) in enumerate(contractions)
            )
            d′ = Diagram(contractions′, d.topology)
            push!(ds′, d′, -coeff)
            contractions′ = FixedVector(
                i ∈ idxs ? make_retarded(e) : e for (i, e) in enumerate(contractions)
            )
            d′ = Diagram(contractions′, d.topology)
            push!(ds′, d′, -coeff)
            continue
        elseif ar_simplification
            contractions′ = FixedVector(
                i ∈ idxs ? make_spectral(e) : e for (i, e) in enumerate(contractions)
            )
            d′ = Diagram(contractions′, d.topology)
            push!(ds′, d′, coeff)
            contractions′ = FixedVector(
                if i == idxs[1]
                    make_retarded(e)
                elseif i == idxs[2]
                    make_advanced(e)
                else
                    e
                end for (i, e) in enumerate(contractions)
            )
            d′ = Diagram(contractions′, d.topology)
            push!(ds′, d′, -coeff)
            continue
        else
            push!(ds′, d′, coeff)
            continue
        end
    end
    filter_nonzero!(ds′)
    return ds′
end
