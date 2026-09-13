"""A canonical independent momentum variable used by exact routing."""
struct MomentumVariable
    index::Int16

    function MomentumVariable(index::Integer)
        index > 0 || throw(ArgumentError("momentum-variable index must be positive"))
        index <= typemax(Int16) || throw(ArgumentError("too many momentum variables"))
        return new(convert(Int16, index))
    end
end

Base.isequal(a::MomentumVariable, b::MomentumVariable) = a.index == b.index
Base.:(==)(a::MomentumVariable, b::MomentumVariable) = isequal(a, b)
Base.hash(x::MomentumVariable, h::UInt) = hash(MomentumVariable, hash(x.index, h))
Base.isless(a::MomentumVariable, b::MomentumVariable) = a.index < b.index

"""Ordered canonical basis for exact linear momenta."""
struct MomentumBasis
    variables::Vector{MomentumVariable}
end

function MomentumBasis(n::Integer)
    n >= 0 || throw(ArgumentError("momentum-basis size must be non-negative"))
    n <= typemax(Int16) || throw(ArgumentError("momentum basis is too large"))
    return MomentumBasis(MomentumVariable[MomentumVariable(i) for i in 1:n])
end

Base.length(basis::MomentumBasis) = length(basis.variables)
Base.getindex(basis::MomentumBasis, i::Int) = basis.variables[i]
Base.iterate(basis::MomentumBasis) = iterate(basis.variables)
Base.iterate(basis::MomentumBasis, state) = iterate(basis.variables, state)
Base.eltype(::Type{MomentumBasis}) = MomentumVariable
Base.IteratorSize(::Type{MomentumBasis}) = Base.HasLength()
Base.isequal(a::MomentumBasis, b::MomentumBasis) = isequal(a.variables, b.variables)
Base.:(==)(a::MomentumBasis, b::MomentumBasis) = isequal(a, b)
Base.hash(x::MomentumBasis, h::UInt) = hash(MomentumBasis, hash(x.variables, h))

const MomentumCoefficient = Rational{Int}

"""Exact linear combination of variables in a `MomentumBasis`."""
struct LinearMomentum
    coefficients::Vector{MomentumCoefficient}
end

function LinearMomentum(coefficients::AbstractVector{<:Integer})
    return LinearMomentum(MomentumCoefficient[c // 1 for c in coefficients])
end

Base.length(momentum::LinearMomentum) = length(momentum.coefficients)
Base.getindex(momentum::LinearMomentum, i::Int) = momentum.coefficients[i]
Base.iszero(momentum::LinearMomentum) = all(iszero, momentum.coefficients)
function Base.zero(momentum::LinearMomentum)
    return LinearMomentum(zeros(MomentumCoefficient, length(momentum)))
end
function Base.isequal(a::LinearMomentum, b::LinearMomentum)
    return isequal(a.coefficients, b.coefficients)
end
Base.:(==)(a::LinearMomentum, b::LinearMomentum) = isequal(a, b)
Base.hash(x::LinearMomentum, h::UInt) = hash(LinearMomentum, hash(x.coefficients, h))

function _check_momentum_dimensions(a::LinearMomentum, b::LinearMomentum)
    length(a) == length(b) || throw(DimensionMismatch("momentum basis sizes differ"))
    return nothing
end

function Base.:+(a::LinearMomentum, b::LinearMomentum)
    _check_momentum_dimensions(a, b)
    out = Vector{MomentumCoefficient}(undef, length(a))
    for i in eachindex(out)
        out[i] = a[i] + b[i]
    end
    return LinearMomentum(out)
end

function Base.:-(a::LinearMomentum)
    out = Vector{MomentumCoefficient}(undef, length(a))
    for i in eachindex(out)
        out[i] = -a[i]
    end
    return LinearMomentum(out)
end

Base.:-(a::LinearMomentum, b::LinearMomentum) = a + (-b)

function Base.:*(coefficient::MomentumCoefficient, momentum::LinearMomentum)
    out = Vector{MomentumCoefficient}(undef, length(momentum))
    for i in eachindex(out)
        out[i] = coefficient * momentum[i]
    end
    return LinearMomentum(out)
end

function Base.:*(coefficient::Integer, momentum::LinearMomentum)
    return (convert(Int, coefficient) // 1) * momentum
end

Base.:*(momentum::LinearMomentum, coefficient::MomentumCoefficient) = coefficient * momentum
Base.:*(momentum::LinearMomentum, coefficient::Integer) = coefficient * momentum

"""Return the `i`th unit momentum in `basis`."""
function basis_momentum(basis::MomentumBasis, i::Int)::LinearMomentum
    checkbounds(basis.variables, i)
    coefficients = zeros(MomentumCoefficient, length(basis))
    coefficients[i] = 1 // 1
    return LinearMomentum(coefficients)
end
