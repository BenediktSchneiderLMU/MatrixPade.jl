# Coefficients c_0, ..., c_N of the matrix power series
#
#   f(x) = exp(i A x + i B x^2) = c_0 + c_1 x + ... ,
#
# built by truncated exponentiation of the matrix polynomial M(x) = i(A x + B x^2).
# Written without any package internals so that it stays an independent check on
# the approximants computed from it.

using LinearAlgebra

function _texp_iab_series(A, B, N)
    T = eltype(A)
    d = size(A, 1)
    zed() = zeros(T, d, d)
    eye() = Matrix{T}(I, d, d)

    M = [zed(), im * A, im * B]

    # truncated product of two coefficient vectors
    function tmul(a, b)
        r = [zed() for _ in 0:N]
        for i in eachindex(a), j in eachindex(b)
            k = (i - 1) + (j - 1)
            k <= N || continue
            r[k+1] += a[i] * b[j]
        end
        return r
    end

    S = [zed() for _ in 0:N]
    S[1] = eye()
    Mk = vcat([eye()], [zed() for _ in 1:N])   # M(x)^k, starting at k = 0
    fact = big(1)
    for k in 1:N
        Mk = tmul(Mk, M)
        fact *= k
        for j in 0:N
            S[j+1] += Mk[j+1] / fact
        end
    end
    return S
end

# The 2x2 pair used throughout the README: A^2 = 2I, B^2 = I, AB = -BA, so
# M(x)^2 = -x^2 (2 + x^2) I and exp(M(x)) has a closed form.
_treadme_A(T) = T[1 im; -im -1]
_treadme_B(T) = T[0 1; 1 0]
