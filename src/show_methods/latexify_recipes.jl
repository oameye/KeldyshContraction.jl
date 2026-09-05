function _postwalk_func(x)
    if x == :𝟙
        return "\\mathbb{1}"
    elseif x == :im
        return :i
    elseif MacroTools.@capture(x, dagger(arg_))
        return "\\bar{$(arg)}"
    elseif MacroTools.@capture(x, power(arg_, reg_))
        return reg > 0 ? "$(arg)^+" : "$(arg)^{-}"
    else
        return x
    end
end

@latexrecipe function f(op::QField)
    cdot --> false
    ex = _to_expression(op)
    ex = MacroTools.postwalk(_postwalk_func, ex)
    return isa(ex, String) ? latexstring(ex) : ex
end

_to_expression(x::Number) = x
function _to_expression(x::Complex)
    iszero(x) && return x
    if iszero(real(x))
        return :($(imag(x))*im)
    elseif iszero(imag(x))
        return real(x)
    else
        return :($(real(x)) + $(imag(x))*im)
    end
end

@latexrecipe function f(op::Union{Diagrams,Diagram,Edge})
    cdot --> false
    return latexify(repr(MIME"text/plain"(), op))
end

function _to_expression(op::Field)
    base = is_barred(op) ? :(dagger($(name(op)))) : :($(name(op)))
    reg = Int(regularisation(op))
    if iszero(reg)
        return base
    elseif reg == 1
        return :(power($base, 1))
    else
        return :(power($base, -1))
    end
end

# Built from the package-native accessors, not `SymbolicUtils.arguments`: the latter
# returns a heterogeneous `Vector{Union{C,Field{S}}}` for interop, which has no place in
# package-owned computation.
function _to_expression(t::QMul)
    fs = fields(t)
    isempty(fs) && return _to_expression(coefficient(t))
    factors = _to_expression.(fs)
    SymbolicUtils._isone(coefficient(t)) ||
        pushfirst!(factors, _to_expression(coefficient(t)))
    return :(*($(factors...)))
end
_to_expression(t::QAdd) = :(+($(_to_expression.(terms(t))...)))
