
###############################
#      BosonicDistribution
###############################

## BosonicDistributionTerm

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
function Base.:*(x::BosonicDistributionTerm, y::BosonicDistributionTerm)
    return BosonicDistributionTerm(vcat(x.momenta, y.momenta))
end

## BosonicDistributions

struct BosonicDistributions
    terms::Dict{BosonicDistributionTerm,ComplexF64}
end
function BosonicDistributions()
    return BosonicDistributions(Dict{BosonicDistributionTerm,ComplexF64}())
end

function Base.push!(
    collection::BosonicDistributions, diagram::BosonicDistributionTerm, prefactor::Number
)
    if haskey(collection.terms, diagram)
        collection.terms[diagram] += prefactor
    else
        collection.terms[diagram] = prefactor
    end
    return collection
end
Base.length(collection::BosonicDistributions) = length(collection.terms)

function Base.:*(term::BosonicDistributions, collection::BosonicDistributions)
    length(term) == 1 || throw(
        ArgumentError("BosonicDistributions can only be multiplied with a single term.")
    )
    out = BosonicDistributions()
    bd, coeff = first(term.terms)
    for (key, val) in collection.terms
        out.terms[bd * key] = coeff * val
    end
    return out
end
function Base.:+(xs::BosonicDistributions, ys::BosonicDistributions)
    out = deepcopy(xs)
    for (key, val) in ys.terms
        push!(out, key, val)
    end
    return out
end

function filter_nonzero!(collection::BosonicDistributions)
    filter!((kv) -> !iszero(kv[2]), collection.terms)
    return collection
end

#################################
#      imaginary_part Im(Σᴿ)
#################################

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
    dict = Dict{keytype(topo),BosonicDistributions}()
    for (t_, keys) in topo
        bda = BosonicDistributions()
        for key in keys
            coeff = ds.diagrams[key]
            bds, coeff′ = imaginary_part(key, coeff)
            push!(bda, bds, coeff′)
        end
        filter_nonzero!(bda)
        dict[t_] = bda
    end
    return dict
end

############################
#          i Σᴷ
############################

"""
⃗a ⃗a = - A A -  ⃗r ⃗r
⃗a ⃖r =  A A -  ⃗r ⃖a
"""
function reduce_to_spectral(ds::Diagrams{E1,E2}) where {E1,E2}
    ds′ = Diagrams{E1,E2}()
    for (d, coeff) in ds.diagrams
        _contractions = contractions(d)
        types = propagator_type.(_contractions)
        directions = direction.(_contractions)

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
                i ∈ idxs ? make_spectral(e) : e for (i, e) in enumerate(_contractions)
            )
            d′ = Diagram(contractions′, d.topology)
            push!(ds′, d′, -coeff)
            contractions′ = FixedVector(
                i ∈ idxs ? make_retarded(e) : e for (i, e) in enumerate(_contractions)
            )
            d′ = Diagram(contractions′, d.topology)
            push!(ds′, d′, -coeff)
            continue
        elseif ar_simplification
            contractions′ = FixedVector(
                i ∈ idxs ? make_spectral(e) : e for (i, e) in enumerate(_contractions)
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
                end for (i, e) in enumerate(_contractions)
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

"""
Convert a `Diagrams` object to a dictionary of `BosonicDistributions`. By substituting:
Gᴷ(k) = 0.5*A

"""

function kelysh_to_distribution(ds::Diagrams)
    topo = topologies(ds)
    dict = Dict{keytype(topo),BosonicDistributions}()
    for (t_, keys) in topo
        bda = BosonicDistributions()
        for key in keys
            coeff = ds.diagrams[key]
            bds, coeff′ = kelysh_to_distribution(key, coeff)
            push!(bda, bds, coeff′)
        end
        filter_nonzero!(bda)
        dict[t_] = bda
    end
    return dict
end

function kelysh_to_distribution(d::Diagram, coeff::ComplexF64=complex(1.0))
    bds = Vector{Momenta}()
    for edge in contractions(d)
        edgetype = propagator_type(edge)
        if is_keldysh(edgetype)
            push!(bds, edge.momenta)
            coeff *= 0.5
        end
    end
    return BosonicDistributionTerm(bds), coeff
end

#################################
#      CollisionIntegral
#################################

struct CollisionIntegral{E}
    terms::Dict{FixedVector{E,Int},BosonicDistributions}
end
function CollisionIntegral(Σ::SelfEnergy{E1,E2}) where {E1,E2}
    # TODO: assert has momenta
    Σk = wigner_transform(Σ)

    tmp = reduce_to_spectral(Σk.keldysh)
    iΣk = kelysh_to_distribution(tmp)

    imΣr = imaginary_part(Σk.retarded)

    Fk = BosonicDistributionTerm([Momenta(0)])
    Fks2 = BosonicDistributions(Dict(Fk => complex(2.0)))

    dict = Dict{FixedVector{E2,Int},BosonicDistributions}()
    for t_ in intersect(keys(iΣk), keys(imΣr))
        dict[t_] = iΣk[t_] + Fks2 * imΣr[t_]
    end
    for t_ in setdiff(keys(imΣr), keys(iΣk))
        dict[t_] = Fks2 * imΣr[t_]
    end
    for t_ in setdiff(keys(iΣk), keys(imΣr))
        dict[t_] = iΣk[t_]
    end
    return CollisionIntegral{E2}(dict)
end
