"""
Exact change of dummy loop-momentum coordinates in one routed momentum basis.

The stored matrix has the convention `p_old = T * p_new` for basis variables, so a routed
linear momentum with row coefficient vector `c` transforms as `c -> c*T`. The physical external
momentum row is fixed. The loop block must be an integer unimodular matrix; loop rows may also
contain an exact integer shift by the external momentum.
"""
struct LoopMomentumTransform
    matrix::Matrix{MomentumCoefficient}
    external_index::Int16

    function LoopMomentumTransform(
        matrix::Matrix{MomentumCoefficient}, external_index::Int16
    )
        size(matrix, 1) == size(matrix, 2) ||
            throw(DimensionMismatch("loop-momentum transform must be square"))
        n = size(matrix, 1)
        ext = Int(external_index)
        1 <= ext <= n || throw(BoundsError(matrix, ext))
        all(value -> denominator(value) == 1, matrix) ||
            throw(ArgumentError("loop-momentum transform entries must be integers"))

        for j in 1:n
            expected = j == ext ? 1 // 1 : 0 // 1
            matrix[ext, j] == expected || throw(
                ArgumentError("external momentum must be fixed by the loop transform")
            )
        end

        loop_indices = Int[i for i in 1:n if i != ext]
        determinant = if isempty(loop_indices)
            one(MomentumCoefficient)
        else
            LinearAlgebra.det(matrix[loop_indices, loop_indices])
        end
        abs(determinant) == 1 ||
            throw(ArgumentError("loop block must be integer unimodular (det = ±1)"))

        return new(copy(matrix), external_index)
    end
end

function LoopMomentumTransform(matrix::AbstractMatrix, external_index::Integer)
    size(matrix, 1) == size(matrix, 2) ||
        throw(DimensionMismatch("loop-momentum transform must be square"))
    n = size(matrix, 1)
    1 <= external_index <= n || throw(BoundsError(matrix, external_index))
    external_index <= typemax(Int16) || throw(
        ArgumentError("loop-momentum basis is too large for the transform index type")
    )

    owned = Matrix{MomentumCoefficient}(undef, n, n)
    for index in eachindex(matrix)
        value = matrix[index]
        (value isa Integer || value isa Rational) || throw(
            ArgumentError("loop-momentum transforms require exact integer/rational data"),
        )
        owned[index] = convert(MomentumCoefficient, value)
    end
    return LoopMomentumTransform(owned, convert(Int16, external_index))
end

loop_transform_matrix(transform::LoopMomentumTransform) = copy(transform.matrix)
external_momentum_index(transform::LoopMomentumTransform) = Int(transform.external_index)

function _external_basis_index(basis::MomentumBasis, external::MomentumVariable)
    matches = Int[i for (i, variable) in enumerate(basis) if variable == external]
    length(matches) == 1 || throw(
        ArgumentError("external momentum must occur exactly once in the momentum basis")
    )
    return only(matches)
end

function _loop_basis_indices(basis::MomentumBasis, external::MomentumVariable)
    external_index = _external_basis_index(basis, external)
    return Int[i for i in eachindex(basis.variables) if i != external_index]
end

"""
Construct an exact loop signed-permutation / external-shift transform.

`permutation[i]` gives the new loop-coordinate slot used by old loop slot `i`. `signs[i]` is
`±1`, and `external_shifts[i]` adds that integer multiple of the fixed external momentum to the
old loop coordinate. The resulting transform is validated through the general unimodular
constructor.
"""
function loop_permutation_transform(
    basis::MomentumBasis,
    external::MomentumVariable,
    permutation::AbstractVector{<:Integer},
    signs::AbstractVector{<:Integer},
    external_shifts::AbstractVector{<:Integer},
)
    external_index = _external_basis_index(basis, external)
    loop_indices = _loop_basis_indices(basis, external)
    nloops = length(loop_indices)
    length(permutation) == nloops ||
        throw(DimensionMismatch("loop permutation has the wrong size"))
    length(signs) == nloops || throw(DimensionMismatch("loop signs have the wrong size"))
    length(external_shifts) == nloops ||
        throw(DimensionMismatch("loop external shifts have the wrong size"))

    sorted_permutation = sort!(Int[convert(Int, value) for value in permutation])
    sorted_permutation == collect(1:nloops) ||
        throw(ArgumentError("loop permutation must contain each loop slot exactly once"))
    all(sign -> sign == 1 || sign == -1, signs) ||
        throw(ArgumentError("loop permutation signs must be ±1"))

    matrix = zeros(MomentumCoefficient, length(basis), length(basis))
    for i in 1:length(basis)
        matrix[i, i] = 1 // 1
    end
    for (slot, old_index) in enumerate(loop_indices)
        for new_index in loop_indices
            matrix[old_index, new_index] = 0 // 1
        end
        matrix[old_index, loop_indices[Int(permutation[slot])]] = Int(signs[slot]) // 1
        matrix[old_index, external_index] = Int(external_shifts[slot]) // 1
    end
    return LoopMomentumTransform(matrix, external_index)
end

function loop_permutation_transform(
    basis::MomentumBasis, external::MomentumVariable, permutation::AbstractVector{<:Integer}
)
    nloops = length(basis) - 1
    return loop_permutation_transform(
        basis, external, permutation, ones(Int8, nloops), zeros(Int, nloops)
    )
end

"""Apply an exact loop-coordinate transform to a routed linear momentum."""
function transform_loop_momenta(
    momentum::LinearMomentum, transform::LoopMomentumTransform
)::LinearMomentum
    n = size(transform.matrix, 1)
    length(momentum) == n ||
        throw(DimensionMismatch("momentum and loop transform use different basis sizes"))
    coefficients = zeros(MomentumCoefficient, n)
    @inbounds for j in 1:n
        value = zero(MomentumCoefficient)
        for i in 1:n
            value += momentum[i] * transform.matrix[i, j]
        end
        coefficients[j] = value
    end
    return LinearMomentum(coefficients)
end

function transform_loop_momenta(
    atom::OccupationAtom{S}, transform::LoopMomentumTransform
) where {S<:Statistics}
    return OccupationAtom{S}(atom.family, transform_loop_momenta(atom.momentum, transform))
end

function transform_loop_momenta(
    monomial::OccupationMonomial{S}, transform::LoopMomentumTransform
) where {S<:Statistics}
    return OccupationMonomial(
        OccupationAtom{S}[transform_loop_momenta(atom, transform) for atom in monomial]
    )
end

function transform_loop_momenta(
    polynomial::OccupationPolynomial{C,S}, transform::LoopMomentumTransform
) where {C<:Number,S<:Statistics}
    terms = Pair{OccupationMonomial{S},C}[
        transform_loop_momenta(monomial, transform) => coefficient for
        (monomial, coefficient) in polynomial
    ]
    return OccupationPolynomial{C,S}(terms)
end

function transform_loop_momenta(
    component::MomentumComponent, transform::LoopMomentumTransform
)
    return MomentumComponent(
        transform_loop_momenta(component.momentum, transform), component.axis
    )
end

function transform_loop_momenta(
    monomial::MomentumMonomial, transform::LoopMomentumTransform
)
    return MomentumMonomial(
        MomentumComponent[
            transform_loop_momenta(component, transform) for component in monomial
        ],
    )
end

function transform_loop_momenta(
    polynomial::MomentumPolynomial{C}, transform::LoopMomentumTransform
) where {C<:Number}
    terms = Pair{MomentumMonomial,C}[
        transform_loop_momenta(monomial, transform) => coefficient for
        (monomial, coefficient) in polynomial
    ]
    return MomentumPolynomial{C}(terms)
end

function transform_loop_momenta(
    atom::DispersionAtom{S}, transform::LoopMomentumTransform
) where {S<:Statistics}
    return DispersionAtom{S}(atom.family, transform_loop_momenta(atom.momentum, transform))
end

function transform_loop_momenta(
    form::EnergyForm{S}, transform::LoopMomentumTransform
) where {S<:Statistics}
    energy_basis_size(form) == size(transform.matrix, 1) ||
        throw(DimensionMismatch("energy form and loop transform use different basis sizes"))
    terms = Pair{DispersionAtom{S},EnergyCoefficient}[
        transform_loop_momenta(atom, transform) => coefficient for
        (atom, coefficient) in form
    ]
    return EnergyForm(energy_basis_size(form), terms)
end

function _transform_frequency_support(
    support::FrequencySupport{S}, transform::LoopMomentumTransform
) where {S<:Statistics}
    shells = EnergyShell{S}[]
    principal_values = PrincipalValueSupport{S}[]
    factor = one(MomentumCoefficient)
    sizehint!(shells, length(support.shells))
    sizehint!(principal_values, length(support.principal_values))

    for shell in support.shells
        transformed, shell_factor = energy_shell(
            transform_loop_momenta(shell.energy, transform)
        )
        push!(shells, transformed)
        factor *= shell_factor
    end
    for principal_value in support.principal_values
        transformed, pv_factor = principal_value_support(
            transform_loop_momenta(principal_value.energy, transform)
        )
        push!(principal_values, transformed)
        factor *= pv_factor
    end
    return FrequencySupport(shells, principal_values), factor
end

"""Canonical regular phase-space sector after the physical loop-coordinate quotient."""
struct CollisionKernelSector{S<:Statistics}
    parameter::ParameterMonomial
    basis::MomentumBasis
    external_momentum::MomentumVariable
    kinematic::MomentumPolynomial{ComplexRationals}
    support::FrequencySupport{S}
end

statistics(::CollisionKernelSector{S}) where {S<:Statistics} = S
parameters(sector::CollisionKernelSector) = sector.parameter
momentum_basis(sector::CollisionKernelSector) = sector.basis
external_wigner_momentum(sector::CollisionKernelSector) = sector.external_momentum
kinematic_factor(sector::CollisionKernelSector) = sector.kinematic
frequency_support(sector::CollisionKernelSector) = sector.support

function Base.isequal(a::CollisionKernelSector{S}, b::CollisionKernelSector{S}) where {S}
    return isequal(a.parameter, b.parameter) &&
           isequal(a.basis, b.basis) &&
           isequal(a.external_momentum, b.external_momentum) &&
           isequal(a.kinematic, b.kinematic) &&
           isequal(a.support, b.support)
end
Base.:(==)(a::CollisionKernelSector, b::CollisionKernelSector) = isequal(a, b)
function Base.hash(sector::CollisionKernelSector, h::UInt)
    h = hash(CollisionKernelSector, h)
    h = hash(sector.parameter, h)
    h = hash(sector.basis, h)
    h = hash(sector.external_momentum, h)
    h = hash(sector.kinematic, h)
    return hash(sector.support, h)
end

function _transform_kernel_sector(
    sector::ReducedCollisionSector{S}, transform::LoopMomentumTransform
) where {S<:Statistics}
    support, factor = _transform_frequency_support(frequency_support(sector), transform)
    transformed = CollisionKernelSector{S}(
        parameters(sector),
        momentum_basis(sector),
        external_wigner_momentum(sector),
        transform_loop_momenta(kinematic_factor(sector), transform),
        support,
    )
    return transformed, factor
end

struct _LoopGraphColor
    kind::UInt8
    name::Symbol
    indices::NTuple{4,Int16}
    axis::Symbol
    a::MomentumCoefficient
    b::MomentumCoefficient
end

const _LOOP_COLOR_NONE = :__loop_none__
const _LOOP_COLOR_NO_INDICES = ntuple(_ -> typemin(Int16), 4)

function Base.isequal(a::_LoopGraphColor, b::_LoopGraphColor)
    return a.kind == b.kind &&
           a.name === b.name &&
           a.indices == b.indices &&
           a.axis === b.axis &&
           a.a == b.a &&
           a.b == b.b
end
Base.:(==)(a::_LoopGraphColor, b::_LoopGraphColor) = isequal(a, b)
function Base.hash(color::_LoopGraphColor, h::UInt)
    h = hash(color.kind, h)
    h = hash(color.name, h)
    h = hash(color.indices, h)
    h = hash(color.axis, h)
    h = hash(color.a, h)
    return hash(color.b, h)
end
function Base.isless(a::_LoopGraphColor, b::_LoopGraphColor)
    a.kind == b.kind || return a.kind < b.kind
    a.name === b.name || return isless(a.name, b.name)
    a.indices == b.indices || return a.indices < b.indices
    a.axis === b.axis || return isless(a.axis, b.axis)
    a.a == b.a || return a.a < b.a
    return a.b < b.b
end

function _loop_graph_color(kind::Integer)
    return _LoopGraphColor(
        UInt8(kind),
        _LOOP_COLOR_NONE,
        _LOOP_COLOR_NO_INDICES,
        _LOOP_COLOR_NONE,
        0 // 1,
        0 // 1,
    )
end

function _loop_family_color(
    kind::Integer,
    family::FieldFamily,
    a::MomentumCoefficient=0 // 1,
    b::MomentumCoefficient=0 // 1,
)
    return _LoopGraphColor(
        UInt8(kind), name(family), slots(field_indices(family)), _LOOP_COLOR_NONE, a, b
    )
end

function _loop_axis_color(kind::Integer, axis::Symbol, a::MomentumCoefficient=0 // 1)
    return _LoopGraphColor(
        UInt8(kind), _LOOP_COLOR_NONE, _LOOP_COLOR_NO_INDICES, axis, a, 0 // 1
    )
end

function _loop_complex_color(kind::Integer, coefficient::ComplexRationals)
    return _LoopGraphColor(
        UInt8(kind),
        _LOOP_COLOR_NONE,
        _LOOP_COLOR_NO_INDICES,
        _LOOP_COLOR_NONE,
        real(coefficient),
        imag(coefficient),
    )
end

mutable struct _LoopCanonicalGraphBuilder
    colors::Vector{_LoopGraphColor}
    edges::Vector{Tuple{Int,Int}}
end

function _LoopCanonicalGraphBuilder()
    return _LoopCanonicalGraphBuilder(_LoopGraphColor[], Tuple{Int,Int}[])
end

function _add_loop_vertex!(builder::_LoopCanonicalGraphBuilder, color::_LoopGraphColor)
    push!(builder.colors, color)
    return length(builder.colors)
end
function _add_loop_edge!(builder::_LoopCanonicalGraphBuilder, source::Int, target::Int)
    push!(builder.edges, (source, target))
    return builder
end

function _add_loop_momentum_incidence!(
    builder::_LoopCanonicalGraphBuilder,
    occurrence::Int,
    momentum::LinearMomentum,
    loop_indices::Vector{Int},
    positive_vertices::Vector{Int},
    negative_vertices::Vector{Int},
    loop_incidences::Vector{Vector{Tuple{Int,MomentumCoefficient}}},
)
    for (slot, basis_index) in enumerate(loop_indices)
        coefficient = momentum[basis_index]
        iszero(coefficient) && continue
        incidence = _add_loop_vertex!(
            builder,
            _LoopGraphColor(
                UInt8(4),
                _LOOP_COLOR_NONE,
                _LOOP_COLOR_NO_INDICES,
                _LOOP_COLOR_NONE,
                abs(coefficient),
                0 // 1,
            ),
        )
        _add_loop_edge!(builder, occurrence, incidence)
        push!(loop_incidences[slot], (occurrence, coefficient))
        orientation = coefficient > 0 ? positive_vertices[slot] : negative_vertices[slot]
        _add_loop_edge!(builder, incidence, orientation)
    end
    return builder
end

function _add_loop_semantics!(
    builder::_LoopCanonicalGraphBuilder,
    root::Int,
    sector::ReducedCollisionSector{S},
    monomial::OccupationMonomial{S},
    external_index::Int,
    loop_indices::Vector{Int},
    positive_vertices::Vector{Int},
    negative_vertices::Vector{Int},
    loop_incidences::Vector{Vector{Tuple{Int,MomentumCoefficient}}},
) where {S<:Statistics}
    for atom in monomial
        momentum = atom.momentum
        vertex = _add_loop_vertex!(
            builder, _loop_family_color(10, atom.family, momentum[external_index])
        )
        _add_loop_edge!(builder, root, vertex)
        _add_loop_momentum_incidence!(
            builder,
            vertex,
            momentum,
            loop_indices,
            positive_vertices,
            negative_vertices,
            loop_incidences,
        )
    end

    for (kinematic_monomial, coefficient) in kinematic_factor(sector)
        term_vertex = _add_loop_vertex!(builder, _loop_complex_color(20, coefficient))
        _add_loop_edge!(builder, root, term_vertex)
        for component in kinematic_monomial
            momentum = component.momentum
            component_vertex = _add_loop_vertex!(
                builder, _loop_axis_color(21, component.axis, momentum[external_index])
            )
            _add_loop_edge!(builder, term_vertex, component_vertex)
            _add_loop_momentum_incidence!(
                builder,
                component_vertex,
                momentum,
                loop_indices,
                positive_vertices,
                negative_vertices,
                loop_incidences,
            )
        end
    end

    support = frequency_support(sector)
    for shell in support.shells
        support_vertex = _add_loop_vertex!(builder, _loop_graph_color(30))
        _add_loop_edge!(builder, root, support_vertex)
        for (atom, coefficient) in shell.energy
            momentum = atom.momentum
            energy_vertex = _add_loop_vertex!(
                builder,
                _loop_family_color(32, atom.family, coefficient, momentum[external_index]),
            )
            _add_loop_edge!(builder, support_vertex, energy_vertex)
            _add_loop_momentum_incidence!(
                builder,
                energy_vertex,
                momentum,
                loop_indices,
                positive_vertices,
                negative_vertices,
                loop_incidences,
            )
        end
    end
    for principal_value in support.principal_values
        support_vertex = _add_loop_vertex!(builder, _loop_graph_color(31))
        _add_loop_edge!(builder, root, support_vertex)
        for (atom, coefficient) in principal_value.energy
            momentum = atom.momentum
            energy_vertex = _add_loop_vertex!(
                builder,
                _loop_family_color(32, atom.family, coefficient, momentum[external_index]),
            )
            _add_loop_edge!(builder, support_vertex, energy_vertex)
            _add_loop_momentum_incidence!(
                builder,
                energy_vertex,
                momentum,
                loop_indices,
                positive_vertices,
                negative_vertices,
                loop_incidences,
            )
        end
    end
    return builder
end

function _loop_graph(builder::_LoopCanonicalGraphBuilder)
    colors = sort!(unique(copy(builder.colors)))
    labels = Int[searchsortedfirst(colors, color) for color in builder.colors]
    graph = NautyGraphs.NautyDiGraph(length(labels); vertex_labels=labels)
    for (source, target) in builder.edges
        Graphs.add_edge!(graph, source, target)
    end
    return graph
end

"""
Return the deterministic signed-permutation gauge of one complete regular integrand atom.

The graph contains the occupation monomial, derivative kinematics, and shell/PV support together.
Each loop variable is represented by an unlabeled orientation pair, so Nauty canonicalization
chooses loop order without enumerating the `2^L L!` signed-permutation group. The sign convention
then makes the earliest canonical semantic occurrence of each loop coordinate positive.
"""
function _canonical_loop_transform(
    sector::ReducedCollisionSector{S}, monomial::OccupationMonomial{S}
) where {S<:Statistics}
    basis = momentum_basis(sector)
    external = external_wigner_momentum(sector)
    external_index = _external_basis_index(basis, external)
    loop_indices = _loop_basis_indices(basis, external)
    nloops = length(loop_indices)

    builder = _LoopCanonicalGraphBuilder()
    root = _add_loop_vertex!(builder, _loop_graph_color(1))
    pair_vertices = Vector{Int}(undef, nloops)
    positive_vertices = Vector{Int}(undef, nloops)
    negative_vertices = Vector{Int}(undef, nloops)
    loop_incidences = [Tuple{Int,MomentumCoefficient}[] for _ in 1:nloops]
    for slot in 1:nloops
        pair = _add_loop_vertex!(builder, _loop_graph_color(2))
        positive = _add_loop_vertex!(builder, _loop_graph_color(3))
        negative = _add_loop_vertex!(builder, _loop_graph_color(3))
        pair_vertices[slot] = pair
        positive_vertices[slot] = positive
        negative_vertices[slot] = negative
        _add_loop_edge!(builder, root, pair)
        _add_loop_edge!(builder, pair, positive)
        _add_loop_edge!(builder, pair, negative)
    end

    _add_loop_semantics!(
        builder,
        root,
        sector,
        monomial,
        external_index,
        loop_indices,
        positive_vertices,
        negative_vertices,
        loop_incidences,
    )

    permutation = NautyGraphs.canonical_permutation(_loop_graph(builder))
    rank = Vector{Int}(undef, length(permutation))
    for (canonical_rank, original_vertex) in enumerate(permutation)
        rank[original_vertex] = canonical_rank
    end

    ordered_slots = sortperm(1:nloops; by=slot -> rank[pair_vertices[slot]])
    loop_permutation = zeros(Int, nloops)
    loop_signs = ones(Int, nloops)
    for (new_slot, old_slot) in enumerate(ordered_slots)
        loop_permutation[old_slot] = new_slot
        incidences = loop_incidences[old_slot]
        isempty(incidences) && continue
        canonical_first = first(sort!(copy(incidences); by=record -> rank[first(record)]))
        loop_signs[old_slot] = last(canonical_first) > 0 ? 1 : -1
    end
    return loop_permutation_transform(
        basis, external, loop_permutation, loop_signs, zeros(Int, nloops)
    )
end

"""
Regular occupation expression after quotienting physically dummy loop-coordinate labels.

Each complete integrand atom is put in one deterministic signed-permutation gauge before exact
coefficient merging. The representation is therefore an actual quotient class, not a Reynolds
average over coordinate spellings.
"""
struct LoopQuotientedExpression{C<:Number,S<:Statistics,O,G,Ctx<:AbstractWignerContext}
    terms::Dict{CollisionKernelSector{S},OccupationPolynomial{C,S}}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::LoopQuotientedExpression{C,S,O}) where {C,S,O} = O
statistics(::LoopQuotientedExpression{C,S}) where {C,S} = S
target_family(expression::LoopQuotientedExpression) = expression.target
parameters(expression::LoopQuotientedExpression) = expression.parameter
gradient_order(::LoopQuotientedExpression{C,S,O,G}) where {C,S,O,G} = Val(G)
wigner_context(expression::LoopQuotientedExpression) = expression.context
loop_quotient_terms(expression::LoopQuotientedExpression) = expression.terms
Base.length(expression::LoopQuotientedExpression) = length(expression.terms)
Base.isempty(expression::LoopQuotientedExpression) = isempty(expression.terms)

"""
Canonical regular quasiparticle phase-space kernel after the physical dummy-loop quotient.

Only regular strict-QP terms enter this representation. Dependent-shell, unresolved causal, and
Trotter branches remain upstream in `ReducedFrequencyCollision`.
"""
struct CollisionKernel{C<:Number,S<:Statistics,O,G,Ctx<:AbstractWignerContext}
    terms::Dict{CollisionKernelSector{S},OccupationPolynomial{C,S}}
    target::FieldFamily{S}
    parameter::ParameterMonomial
    context::Ctx
end

order(::CollisionKernel{C,S,O}) where {C,S,O} = O
statistics(::CollisionKernel{C,S}) where {C,S} = S
target_family(kernel::CollisionKernel) = kernel.target
parameters(kernel::CollisionKernel) = kernel.parameter
gradient_order(::CollisionKernel{C,S,O,G}) where {C,S,O,G} = Val(G)
wigner_context(kernel::CollisionKernel) = kernel.context
collision_kernel_terms(kernel::CollisionKernel) = kernel.terms
Base.length(kernel::CollisionKernel) = length(kernel.terms)
Base.isempty(kernel::CollisionKernel) = isempty(kernel.terms)

function _push_kernel_polynomial!(
    terms::Dict{CollisionKernelSector{S},OccupationPolynomial{C,S}},
    sector::CollisionKernelSector{S},
    polynomial::OccupationPolynomial{C,S},
) where {C<:Number,S<:Statistics}
    iszero(polynomial) && return terms
    if haskey(terms, sector)
        combined = terms[sector] + polynomial
        if iszero(combined)
            delete!(terms, sector)
        else
            terms[sector] = combined
        end
    else
        terms[sector] = polynomial
    end
    return terms
end

function _kinematic_atom_sector(
    sector::ReducedCollisionSector{S}, monomial::MomentumMonomial
) where {S<:Statistics}
    kinematic = MomentumPolynomial(monomial, one(ComplexRationals))
    return ReducedCollisionSector{S}(
        parameters(sector),
        momentum_basis(sector),
        external_wigner_momentum(sector),
        kinematic,
        frequency_support(sector),
    )
end

"""
    quotient_loop_momenta(expression)

Canonicalize each regular integrand atom under signed permutations of dummy loop coordinates and
merge equivalent classes exactly. Before choosing the loop gauge, each derivative kinematic
polynomial is distributed into unit-coefficient momentum monomials. Scalar kinematic coefficients
are folded into the occupation coefficient, so factored and expanded phase-space expressions have
the same canonical representation. Occupation factors, derivative kinematics, and shell/PV
support are transformed coherently. The canonical gauge is obtained by graph canonicalization
rather than an explicit factorial search.

General exact unimodular transforms and external-momentum shifts remain available through
`LoopMomentumTransform`; automatic quotienting currently uses the universally safe signed-
permutation subgroup.
"""
function quotient_loop_momenta(
    expression::OccupationReducedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:Statistics,O,G,Ctx<:AbstractWignerContext}
    D = promote_type(C, ComplexRationals, Rational{Int})
    out = Dict{CollisionKernelSector{S},OccupationPolynomial{D,S}}()

    for (sector, polynomial) in occupation_reduced_terms(expression)
        for (occupation_monomial, occupation_coefficient) in polynomial
            for (kinematic_monomial, kinematic_coefficient) in kinematic_factor(sector)
                atom_sector = _kinematic_atom_sector(sector, kinematic_monomial)
                transform = _canonical_loop_transform(atom_sector, occupation_monomial)
                transformed_sector, support_factor = _transform_kernel_sector(
                    atom_sector, transform
                )
                transformed_monomial = transform_loop_momenta(
                    occupation_monomial, transform
                )
                transformed_coefficient =
                    convert(D, occupation_coefficient) *
                    convert(D, kinematic_coefficient) *
                    convert(D, support_factor)
                contribution = Pair{OccupationMonomial{S},D}[transformed_monomial => transformed_coefficient]
                _push_kernel_polynomial!(
                    out, transformed_sector, OccupationPolynomial{D,S}(contribution)
                )
            end
        end
    end

    return LoopQuotientedExpression{D,S,O,G,Ctx}(
        out, target_family(expression), parameters(expression), wigner_context(expression)
    )
end

"""Lower a loop-quotiented regular expression to the final symbolic collision-kernel IR."""
function collision_kernel(
    expression::LoopQuotientedExpression{C,S,O,G,Ctx}
) where {C<:Number,S<:Statistics,O,G,Ctx<:AbstractWignerContext}
    return CollisionKernel{C,S,O,G,Ctx}(
        expression.terms,
        target_family(expression),
        parameters(expression),
        wigner_context(expression),
    )
end

"""Quotient dummy loop coordinates and construct the regular strict-QP collision kernel."""
function collision_kernel(expression::OccupationReducedExpression)
    return collision_kernel(quotient_loop_momenta(expression))
end
