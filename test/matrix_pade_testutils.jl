# Small polynomial helpers for the matrix Pade tests, written without using
# any package internals so that the checks stay independent of the code they
# are testing.

# Degree of a coefficient vector; -1 stands in for the zero polynomial.
_tdeg(a) = something(findlast(!iszero, a), 0) - 1

# Drop trailing zeros, so comparisons do not depend on how a polynomial
# happens to be padded.
_ttrim(a) = a[1:_tdeg(a)+1]

# dct P = min_lam (n_lam + 1 - deg P_lam), skipping the zero components
# (Beckermann and Labahn 1994, Definition 3.1).
function _tdefect(row, n)
    ds = [n[lam] + 1 - _tdeg(row[lam]) for lam in eachindex(row) if _tdeg(row[lam]) >= 0]
    return isempty(ds) ? typemax(Int) : minimum(ds)
end

# The s-residual P_1(z^s) f_1(z) + ... + P_m(z^s) f_m(z), truncated to
# degree sigma-1. All zero means ord P >= sigma.
function _tresidual(row, F, s, sigma)
    acc = zeros(eltype(F[1]), sigma)
    for lam in eachindex(row), i in eachindex(row[lam])
        a = row[lam][i]
        iszero(a) && continue
        for k in 1:sigma
            idx = k + s * (i - 1)
            idx <= sigma || break
            acc[idx] += a * F[lam][k]
        end
    end
    return acc
end

# Coefficient of z^j in A(z) Q(z) - P(z) (side = :right) or in
# Q(z) A(z) - P(z) (side = :left), by plain truncated polynomial arithmetic.
function _torder_coeff(A, Pc, Qc, j, side)
    acc = zero(side === :right ? A[1] * Qc[1] : Qc[1] * A[1])
    for k in 0:min(j, length(Qc) - 1)
        acc = acc + (side === :right ? A[j-k+1] * Qc[k+1] : Qc[k+1] * A[j-k+1])
    end
    j < length(Pc) && (acc = acc - Pc[j+1])
    return acc
end
