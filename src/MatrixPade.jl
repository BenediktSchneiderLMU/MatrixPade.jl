module MatrixPade

using LinearAlgebra

export matrix_pade_eval

"""
    matrix_pade_eval(coeffs, x)

Placeholder stub. Will evaluate the matrix-Padé approximation defined by
`coeffs` (a vector of matrix-valued coefficients) at `x`. Currently just
returns the leading coefficient unchanged, pending the real implementation.
"""
function matrix_pade_eval(coeffs::AbstractVector{<:AbstractMatrix}, x)
    isempty(coeffs) && throw(ArgumentError("coeffs must be non-empty"))
    return first(coeffs)
end

end # module MatrixPade
