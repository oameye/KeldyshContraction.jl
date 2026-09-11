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

struct BosonicDistributions{C<:Number}
    terms::Dict{BosonicDistributionTerm,C}
end
function BosonicDistributions{C}() where {C<:Number}
    return BosonicDistributions{C}(Dict{BosonicDistributionTerm,C}())
end
BosonicDistributions() = BosonicDistributions{ComplexRationals}()

function Base.push!(
    collection::BosonicDistributions{C}, diagram::BosonicDistributionTerm, prefactor::Number
) where {C<:Number}
    value = _simplify(convert(C, prefactor))
    if haskey(collection.terms, diagram)
        collection.terms[diagram] = _simplify(collection.terms[diagram] + value)
    else
        collection.terms[diagram] = value
    end
    return collection
end
Base.length(collection::BosonicDistributions) = length(collection.terms)
Base.isempty(collection::BosonicDistributions) = isempty(collection.terms)

function Base.:*(
    term::BosonicDistributions{C1}, collection::BosonicDistributions{C2}
) where {C1<:Number,C2<:Number}
    length(term) == 1 || throw(
        ArgumentError("BosonicDistributions can only be multiplied with a single term.")
    )
    D = promote_type(C1, C2)
    out = BosonicDistributions{D}()
    bd, coeff = first(term.terms)
    coeff′ = convert(D, coeff)
    for (key, val) in collection.terms
        out.terms[bd * key] = _simplify(coeff′ * convert(D, val))
    end
    return out
end
function Base.:*(x::P, collection::BosonicDistributions{C}) where {P<:Number,C<:Number}
    D = promote_type(P, C)
    out = BosonicDistributions{D}()
    x′ = convert(D, x)
    for (key, val) in collection.terms
        out.terms[key] = _simplify(x′ * convert(D, val))
    end
    return out
end
function Base.:+(
    xs::BosonicDistributions{C1}, ys::BosonicDistributions{C2}
) where {C1<:Number,C2<:Number}
    D = promote_type(C1, C2)
    out = BosonicDistributions{D}()
    for (key, val) in xs.terms
        push!(out, key, val)
    end
    for (key, val) in ys.terms
        push!(out, key, val)
    end
    return out
end

function filter_nonzero!(collection::BosonicDistributions)
    filter!((kv) -> !iszero(kv[2]), collection.terms)
    return collection
end

function _real_distribution_coefficient_type(::Type{C}) where {C<:Number}
    return promote_type(C, Rational{Int64})
end
function _distribution_coefficient_type(::Type{C}) where {C<:Number}
    return promote_type(C, ComplexRationals)
end

#################################
#      imaginary_part Im(Σᴿ)
#################################
"""
Im(Gᴿ) = -0.5 * A
Im(Gᴬ) = 0.5 * A
Im(Gᴷ) = 0.5 * F * A
"""
function imaginary_part(d::Diagram, coeff::C) where {C<:Number}
    D = _real_distribution_coefficient_type(C)
    value = convert(D, coeff)
    half = convert(D, 1 // 2)
    bds = Vector{Momenta}()
    for edge in contractions(d)
        edgetype = propagator_type(edge)
        if is_keldysh(edgetype)
            push!(bds, edge.momenta)
        end
        value *= is_advanced(edgetype) ? half : -half
    end
    return BosonicDistributionTerm(bds), _simplify(value)
end
imaginary_part(d::Diagram) = imaginary_part(d, one(ComplexRationals))

function imaginary_part(ds::Diagrams{C,S,E1,E2}) where {C<:Number,S<:Statistics,E1,E2}
    D = _real_distribution_coefficient_type(C)
    topo = topologies(ds)
    dict = Dict{FixedVector{E2,Int},BosonicDistributions{D}}()
    for (t_, keys) in topo
        bda = BosonicDistributions{D}()
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
function reduce_to_spectral(ds::Diagrams{C,S,E1,E2}) where {C<:Number,S<:Statistics,E1,E2}
    ds′ = Diagrams{C,S,E1,E2}()
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
            contractions′ = FixedVector{E1,Edge{S}}(
                i ∈ idxs ? make_spectral(e) : e for (i, e) in enumerate(_contractions)
            )
            d′ = Diagram(contractions′, d.topology)
            push!(ds′, d′, -coeff)
            contractions′ = FixedVector{E1,Edge{S}}(
                i ∈ idxs ? make_retarded(e) : e for (i, e) in enumerate(_contractions)
            )
            d′ = Diagram(contractions′, d.topology)
            push!(ds′, d′, -coeff)
            continue
        elseif ar_simplification
            contractions′ = FixedVector{E1,Edge{S}}(
                i ∈ idxs ? make_spectral(e) : e for (i, e) in enumerate(_contractions)
            )
            d′ = Diagram(contractions′, d.topology)
            push!(ds′, d′, coeff)
            contractions′ = FixedVector{E1,Edge{S}}(
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
        end
    end
    filter_nonzero!(ds′)
    return ds′
end

"""
Convert a `Diagrams` object to a dictionary of `BosonicDistributions`. By substituting:
Gᴷ(k) = im*0.5*F(k)A(k)

"""
function kelysh_to_distribution(
    ds::Diagrams{C,S,E1,E2}
) where {C<:Number,S<:Statistics,E1,E2}
    D = _distribution_coefficient_type(C)
    topo = topologies(ds)
    dict = Dict{FixedVector{E2,Int},BosonicDistributions{D}}()
    for (t_, keys) in topo
        bda = BosonicDistributions{D}()
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

function kelysh_to_distribution(d::Diagram, coeff::C) where {C<:Number}
    D = _distribution_coefficient_type(C)
    value = convert(D, coeff)
    factor = convert(D, (1 // 2) * im)
    bds = Vector{Momenta}()
    for edge in contractions(d)
        edgetype = propagator_type(edge)
        if is_keldysh(edgetype)
            push!(bds, edge.momenta)
            value *= factor
        end
    end
    return BosonicDistributionTerm(bds), _simplify(value)
end
kelysh_to_distribution(d::Diagram) = kelysh_to_distribution(d, one(ComplexRationals))

#################################
#      CollisionIntegral
#################################

struct CollisionIntegral{C<:Number,E}
    terms::Dict{FixedVector{E,Int},BosonicDistributions{C}}
end
function CollisionIntegral(Σ::SelfEnergy{C,S,O,E1,E2}) where {C<:Number,S,O,E1,E2}
    D = _distribution_coefficient_type(C)
    Σk = _legacy_momentum_transform(Σ)

    tmp = reduce_to_spectral(Σk.keldysh)
    ΣkF = kelysh_to_distribution(tmp)

    imΣr = imaginary_part(Σk.retarded)

    Fk = BosonicDistributionTerm([Momenta(0)])
    Fks2 = BosonicDistributions{D}(Dict(Fk => convert(D, 2)))

    dict = Dict{FixedVector{E2,Int},BosonicDistributions{D}}()
    for t_ in intersect(keys(ΣkF), keys(imΣr))
        dict[t_] = im * ΣkF[t_] + Fks2 * imΣr[t_]
    end
    for t_ in setdiff(keys(imΣr), keys(ΣkF))
        dict[t_] = Fks2 * imΣr[t_]
    end
    for t_ in setdiff(keys(ΣkF), keys(imΣr))
        dict[t_] = im * ΣkF[t_]
    end
    return CollisionIntegral{D,E2}(dict)
end
