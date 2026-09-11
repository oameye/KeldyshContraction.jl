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
