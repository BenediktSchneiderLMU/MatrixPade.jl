# Matrix Padé-type approximant (MPTA), scalar denominator / matrix
# numerator case: Gu (2004), Theorem 4.3 (which subsumes Theorem 4.2 as the
# special case m = n-1). 

_exactify(coeffs::AbstractVector{<:AbstractMatrix{<:Integer}}) =
    [Rational{BigInt}.(c) for c in coeffs]
_exactify(coeffs::AbstractVector{<:AbstractMatrix}) = coeffs

"""
    mpta_coeffs(coeffs, m, n) -> (qcoeffs, Pcoeffs)

Compute the matrix Padé-type approximant `(m/n)_f` of Gu (2004) for the
matrix power series `f(z) = coeffs[1] + coeffs[2] z + ...`, with a scalar
denominator `q(z)` and matrix-valued numerator `P(z)` satisfying
`q(z) f(z) - P(z) = O(z^{m+1})`.

Returns `qcoeffs = [q0, ..., qn]` (so `q(z) = q0 + q1 z + ... + qn z^n`)
and `Pcoeffs = [P0, ..., Pm]` (so `P(z) = P0 + P1 z + ... + Pm z^m`).

Requires `n >= 1`, `m >= n-1`, and at least `m+n+1` coefficients. If the
coefficient matrices have an `Integer` element type they are promoted to
`Rational{BigInt}` so the computation is exact.
"""
function mpta_coeffs(coeffs::AbstractVector{<:AbstractMatrix}, m::Integer, n::Integer)
    n >= 1 || throw(ArgumentError("n must be >= 1, got n=$n"))
    m >= n - 1 || throw(ArgumentError("m must be >= n-1, got m=$m, n=$n"))
    length(coeffs) >= m + n + 1 ||
        throw(ArgumentError("need at least m+n+1 = $(m + n + 1) coefficients, got $(length(coeffs))"))

    c = _exactify(coeffs)
    p = m - n + 1 # 0-based index of the first series coefficient used in H

    # c[i+1] holds the 0-based coefficient c_i.
    H = [matdot(c[p+r+1], c[p+r+k+1]) for r in 0:n-1, k in 0:n-1]
    xi = [matdot(c[p+r+1], c[p+r+n+1]) for r in 0:n-1]
    A = hcat(H, xi) # n x (n+1), purely scalar

    S = Vector{eltype(H)}(undef, n + 1)
    for cc in 0:n
        cols = setdiff(0:n, cc) .+ 1
        S[cc+1] = (-1)^(n + cc) * det(A[:, cols])
    end

    qcoeffs = [S[n-j+1] for j in 0:n]

    all(iszero, qcoeffs) &&
        throw(ArgumentError("the Padé-type denominator is identically zero: the Hankel " *
                            "matrix of matdot products is singular for (m, n) = ($m, $n). " *
                            "Gu's scalar product is bilinear and does not conjugate, so it " *
                            "can be isotropic on a complex series; try a different (m, n), " *
                            "or use `matrix_pade` for a genuine matrix denominator."))

    Pcoeffs = Vector{typeof(c[1])}(undef, m + 1)
    for j in 0:m
        acc = zero(c[1])
        for k in 0:min(j, n)
            acc = acc + qcoeffs[k+1] * c[j-k+1]
        end
        Pcoeffs[j+1] = acc
    end

    return qcoeffs, Pcoeffs
end

"""
    PadeApproximant(q, P)

A Padé-style approximant with a scalar-coefficient denominator `q` and a
matrix-coefficient numerator `P`. Calling `pa(z)` evaluates `P(z) / q(z)`.
"""
struct PadeApproximant{Tq,TP<:AbstractMatrix}
    q::Vector{Tq}
    P::Vector{TP}
end

function _horner_scalar(coefs::AbstractVector, z)
    acc = coefs[end]
    for i in length(coefs)-1:-1:1
        acc = acc * z + coefs[i]
    end
    return acc
end

function _horner_matrix(coefs::AbstractVector{<:AbstractMatrix}, z)
    acc = coefs[end]
    for i in length(coefs)-1:-1:1
        acc = acc * z + coefs[i]
    end
    return acc
end

function (pa::PadeApproximant)(z)
    return _horner_matrix(pa.P, z) ./ _horner_scalar(pa.q, z)
end

"""
    mpta(coeffs, m, n) -> PadeApproximant
    mpta(coeffs, m, n, z) -> matrix

Build the `(m/n)` matrix Padé-type approximant (see [`mpta_coeffs`](@ref))
as a [`PadeApproximant`](@ref). If `z` is given, evaluate it immediately
and return the resulting matrix instead.
"""
function mpta(coeffs::AbstractVector{<:AbstractMatrix}, m::Integer, n::Integer)
    qcoeffs, Pcoeffs = mpta_coeffs(coeffs, m, n)
    return PadeApproximant(qcoeffs, Pcoeffs)
end

mpta(coeffs::AbstractVector{<:AbstractMatrix}, m::Integer, n::Integer, z) =
    mpta(coeffs, m, n)(z)

