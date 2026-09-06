"""One factor in a perturbation-parameter monomial."""
struct ParameterPower
    name::Symbol
    exponent::Int

    function ParameterPower(name::Symbol, exponent::Integer)
        exponent > 0 || throw(ArgumentError("parameter exponents must be positive"))
        return new(name, Int(exponent))
    end
end

function Base.isequal(a::ParameterPower, b::ParameterPower)
    return a.name === b.name && a.exponent == b.exponent
end
Base.:(==)(a::ParameterPower, b::ParameterPower) = isequal(a, b)
Base.hash(p::ParameterPower, h::UInt) = hash(p.name, hash(p.exponent, h))
Base.isless(a::ParameterPower, b::ParameterPower) = isless(a.name, b.name)

"""
Canonical commutative monomial in perturbation parameters.

Factors are sorted by name and equal factors are combined, so multiplication order does
not affect equality or hashing. This is the package-native perturbative bookkeeping type;
SymbolicUtils expressions are converted to it at API boundaries.
"""
struct ParameterMonomial
    powers::Vector{ParameterPower}

    function ParameterMonomial(powers::Vector{ParameterPower})
        isempty(powers) && return new(ParameterPower[])
        sorted = sort(copy(powers))
        combined = ParameterPower[]
        current = first(sorted)
        for power in Iterators.drop(sorted, 1)
            if power.name === current.name
                current = ParameterPower(current.name, current.exponent + power.exponent)
            else
                push!(combined, current)
                current = power
            end
        end
        push!(combined, current)
        return new(combined)
    end
end

ParameterMonomial(name::Symbol) = ParameterMonomial(ParameterPower[ParameterPower(name, 1)])

function Base.isequal(a::ParameterMonomial, b::ParameterMonomial)
    return isequal(a.powers, b.powers)
end
Base.:(==)(a::ParameterMonomial, b::ParameterMonomial) = isequal(a, b)
Base.hash(p::ParameterMonomial, h::UInt) = hash(ParameterMonomial, hash(p.powers, h))
Base.one(::Type{ParameterMonomial}) = ParameterMonomial(ParameterPower[])
Base.one(::ParameterMonomial) = one(ParameterMonomial)
Base.isone(p::ParameterMonomial) = isempty(p.powers)

function Base.:*(a::ParameterMonomial, b::ParameterMonomial)
    return ParameterMonomial(vcat(a.powers, b.powers))
end
function Base.:^(p::ParameterMonomial, exponent::Integer)
    exponent >= 0 || throw(DomainError(exponent, "parameter powers must be non-negative"))
    iszero(exponent) && return one(ParameterMonomial)
    return ParameterMonomial(
        ParameterPower[
            ParameterPower(power.name, power.exponent * exponent) for power in p.powers
        ],
    )
end

function Base.show(io::IO, p::ParameterMonomial)
    isempty(p.powers) && return write(io, "1")
    for (i, power) in enumerate(p.powers)
        i > 1 && write(io, "*")
        write(io, string(power.name))
        power.exponent == 1 || write(io, "^", string(power.exponent))
    end
    return nothing
end

"""
    parameter_monomial(x)

Normalize a perturbation parameter to a canonical `ParameterMonomial`.

Symbols, products, and non-negative integer powers from SymbolicUtils are accepted as an
interop boundary. Numeric input is accepted only for the multiplicative identity `1`.
"""
parameter_monomial(p::ParameterMonomial) = p
parameter_monomial(name::Symbol) = ParameterMonomial(name)
function parameter_monomial(x::Number)
    isone(x) ||
        throw(ArgumentError("numeric perturbation parameters are only supported for 1"))
    return one(ParameterMonomial)
end

function parameter_monomial(x::CSym)
    if !SymbolicUtils.iscall(x)
        return ParameterMonomial(Symbol(string(x)))
    end

    operation = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if operation === (*)
        result = one(ParameterMonomial)
        for arg in args
            result *= parameter_monomial(arg)
        end
        return result
    elseif operation === (^)
        length(args) == 2 || throw(ArgumentError("invalid parameter power"))
        exponent = args[2]
        exponent isa Integer || throw(ArgumentError("parameter exponents must be integers"))
        return parameter_monomial(args[1])^exponent
    end
    return throw(ArgumentError("perturbation parameters must be symbols, products, or powers"))
end
