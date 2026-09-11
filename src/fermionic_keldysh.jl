"""Fermionic field statistics."""
struct Fermion <: Statistics end

"""First fermionic Larkin-Ovchinnikov component."""
const One = KeldyshIndex.First

"""Second fermionic Larkin-Ovchinnikov component."""
const Two = KeldyshIndex.Second

is_one(f::Field{Fermion}) = keldysh_index(f) === One
is_two(f::Field{Fermion}) = keldysh_index(f) === Two

"""Exchange of two fermionic generators contributes a minus sign."""
exchange_sign(::Type{Fermion}) = Int8(-1)

@inline function has_duplicate_fermion(args::Vector{Field{Fermion}})
    @inbounds for i in 2:length(args)
        isequal(args[i - 1], args[i]) && return true
    end
    return false
end

# Fermionic products use the shared QMul storage, but canonicalization must additionally
# enforce Grassmann nilpotency. The result remains QMul{C,Fermion} for all runtime values.
@inline function canonical_qmul(arg_c::C, args_nc::Vector{Field{Fermion}}) where {C<:Number}
    iszero(arg_c) && return QMul{C,Fermion}(zero(C), Field{Fermion}[], Val(:presorted))

    sign = canonicalize_fields!(args_nc)
    has_duplicate_fermion(args_nc) &&
        return QMul{C,Fermion}(zero(C), Field{Fermion}[], Val(:presorted))

    coeff = sign == 1 ? arg_c : -arg_c
    return QMul{C,Fermion}(coeff, args_nc, Val(:presorted))
end

"""Permutation parity used by fermionic Wick pairings."""
function pairing_sign(::Type{Fermion}, permutation)::Int8
    sign = Int8(1)
    @inbounds for i in 1:(length(permutation) - 1)
        for j in (i + 1):length(permutation)
            permutation[i] > permutation[j] && (sign = -sign)
        end
    end
    return sign
end

function contraction_compatible(::Type{Fermion}, out::Field{Fermion}, in::Field{Fermion})
    return same_field_family(out, in)
end

@inline function fermionic_structural_zero(out::Field{Fermion}, in::Field{Fermion})
    return is_two(out) && is_one(in)
end

function contraction_filter(v::Contraction{Fermion})
    out, in = v
    contraction_compatible(Fermion, out, in) || return false
    balanced_orientation(v) || return false
    is_physical_propagator(v) || return false
    return !fermionic_structural_zero(out, in)
end

function propagator_checks(
    ::Type{Fermion}, out::Field{Fermion}, in::Field{Fermion}
)::Nothing
    @assert is_barred(in) "The incoming fermionic field must be barred"
    @assert is_unbarred(out) "The outgoing fermionic field must be unbarred"
    @assert contraction_compatible(Fermion, out, in) "Contracted fermionic fields must belong to the same field family"

    v = Contraction(out, in)
    ps = position.(v)
    @assert !is_in(first(ps)) "The outgoing field cannot be at In()"
    @assert !is_out(last(ps)) "The incoming field cannot be at Out()"
    @assert !(has_in(ps) && has_out(ps)) "Cannot contract In() directly with Out()"
    @assert !fermionic_structural_zero(out, in) "The fermionic (Two, One) propagator is zero"
    return nothing
end
function propagator_checks(out::Field{Fermion}, in::Field{Fermion})
    return propagator_checks(Fermion, out, in)
end

"""Determine the fermionic LO propagator type."""
function propagator_type(
    ::Type{Fermion}, out::Field{Fermion}, in::Field{Fermion}
)::PropagatorType.T
    if is_one(out)
        return is_one(in) ? PropagatorType.Retarded : PropagatorType.Keldysh
    elseif is_two(in)
        return PropagatorType.Advanced
    else
        throw(ArgumentError("the fermionic (Two, One) propagator is structurally zero"))
    end
end
propagator_type(out::Field{Fermion}, in::Field{Fermion}) = propagator_type(Fermion, out, in)

@inline fermionic_dual_index(index::KeldyshIndex.T) = index === One ? Two : One

function fermionic_adjoint_field(field::Field{Fermion}, position::Position)
    toggled = bar(field)(position)
    return reconstruct(toggled; keldysh=fermionic_dual_index(keldysh_index(field)))
end

function Base.adjoint(c::Contraction{Fermion})
    return Contraction(
        fermionic_adjoint_field(c.in, position(c.out)),
        fermionic_adjoint_field(c.out, position(c.in)),
    )
end
function Base.adjoint(e::Edge{Fermion})
    c = adjoint(Contraction(e.out, e.in))
    return Edge{Fermion}(c.out, c.in, propagator_type(Fermion, c.out, c.in), e.momenta)
end

function structural_zero(
    ::Type{Fermion}, ::Type{Diagrams{C,Fermion,E1,E2}}
) where {C<:Number,E1,E2}
    return Diagrams{C,Fermion,E1,E2}()
end

"""Return the fermionic triangular LO propagator matrix `[[R,K],[0,A]]`."""
function matrix(G::DressedPropagator{C,Fermion,O,E1,E2}) where {C<:Number,O,E1,E2}
    return matrix(Fermion, G)
end
function matrix(
    ::Type{Fermion}, G::DressedPropagator{C,Fermion,O,E1,E2}
) where {C<:Number,O,E1,E2}
    D = Diagrams{C,Fermion,E1,E2}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = G.retarded
    result[1, 2] = G.keldysh
    result[2, 1] = structural_zero(Fermion, D)
    result[2, 2] = G.advanced
    return result
end

function propagator_fields(L::InteractionLagrangian{C,Fermion}, ::Nothing) where {C<:Number}
    family = target_family(L)
    return family[One], family[Two]
end
function propagator_fields(
    L::InteractionLagrangian{C,Fermion}, target::FieldFamily{Fermion}
) where {C<:Number}
    family = target_family(L, target)
    return family[One], family[Two]
end

function DressedPropagator(
    L::InteractionLagrangian{C,Fermion},
    ::Val{O},
    ::Val{E};
    target=nothing,
    simplify=true,
    _set_reg_to_zero=true,
    kwargs...,
) where {C<:Number,O,E}
    @assert number_of_propagators(L) * O + 1 == E "The supplied Val{edges} must equal the interaction's propagator count times Val{order}, plus the external propagator"
    onefield, twofield = propagator_fields(L, target)

    keldysh = wick_contraction(
        onefield(Out()) * bar(twofield)(In()),
        L,
        Val(O),
        Val(E);
        simplify,
        _set_reg_to_zero,
        kwargs...,
    )
    retarded = wick_contraction(
        onefield(Out()) * bar(onefield)(In()),
        L,
        Val(O),
        Val(E);
        simplify,
        _set_reg_to_zero,
        kwargs...,
    )
    advanced = wick_contraction(
        twofield(Out()) * bar(twofield)(In()),
        L,
        Val(O),
        Val(E);
        simplify,
        _set_reg_to_zero,
        kwargs...,
    )

    for component in (keldysh, retarded, advanced)
        filter_nonzero!(component)
    end

    return DressedPropagator(keldysh, retarded, advanced, Val(O), parameters(L)^O)
end

function self_energy_type(::Type{Fermion}, dict::SmallCollections.SmallDict)
    if is_keldysh(dict[:out]) && is_advanced(dict[:in])
        return PropagatorType.Advanced
    elseif is_retarded(dict[:out]) && is_keldysh(dict[:in])
        return PropagatorType.Retarded
    elseif is_retarded(dict[:out]) && is_advanced(dict[:in])
        return PropagatorType.Keldysh
    else
        error("Fermionic LO self-energy external propagator combination is zero")
    end
end

function _self_energy(G::DressedPropagator{C,Fermion,O,E1,E2}) where {C<:Number,O,E1,E2}
    SE = E1 - 2
    ST = max_edges(O)
    D = Diagrams{C,Fermion,SE,ST}
    self_energy = SmallCollections.SmallDict{3,PropagatorType.T,D}((
        PropagatorType.Advanced => D(),
        PropagatorType.Retarded => D(),
        PropagatorType.Keldysh => D(),
    ))
    construct_self_energy!(self_energy, G.keldysh)

    _simplify_prefactors!(self_energy[PropagatorType.Keldysh])
    _simplify_prefactors!(self_energy[PropagatorType.Retarded])
    _simplify_prefactors!(self_energy[PropagatorType.Advanced])
    return SelfEnergy{C,Fermion,O,SE,ST}(
        self_energy[PropagatorType.Keldysh],
        self_energy[PropagatorType.Retarded],
        self_energy[PropagatorType.Advanced],
        G.parameter,
    )
end

"""Return the fermionic triangular LO self-energy matrix `[[R,K],[0,A]]`."""
function matrix(Σ::SelfEnergy{C,Fermion,O,E1,E2}) where {C<:Number,O,E1,E2}
    return matrix(Fermion, Σ)
end
function matrix(::Type{Fermion}, Σ::SelfEnergy{C,Fermion,O,E1,E2}) where {C<:Number,O,E1,E2}
    D = Diagrams{C,Fermion,E1,E2}
    result = Matrix{D}(undef, 2, 2)
    result[1, 1] = Σ.retarded
    result[1, 2] = Σ.keldysh
    result[2, 1] = structural_zero(Fermion, D)
    result[2, 2] = Σ.advanced
    return result
end
