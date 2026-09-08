# Fast power Hermite Padé solver (FPHPS): Beckermann and Labahn (1994), §3,
# p. 810. Computes a sigma-basis of the power Hermite Padé approximants
# (PHPAs) of type (n, sigma, s) for a vector F of scalar power series.

_fphps_eltype(::Type{<:Integer}) = Rational{BigInt}
_fphps_eltype(::Type{T}) where {T} = T

# a .-= f * b, growing `a` with zeros first if `b` is the longer polynomial.
function _polysub!(a::Vector{T}, b::AbstractVector, f) where {T}
    if length(b) > length(a)
        append!(a, zeros(T, length(b) - length(a)))
    end
    @inbounds for k in eachindex(b)
        a[k] -= f * b[k]
    end
    return a
end

"""
    fphps(F, sigma, n, s) -> (P, d, pivots)

Compute a `sigma`-basis of power Hermite Padé approximants (PHPAs) of type
`(n, sigma, s)` for the scalar power series vector `F = (f_1, ..., f_m)`,
following the FPHPS algorithm of Beckermann and Labahn (1994), §3.

A PHPA is a row of polynomials `P = (P_1, ..., P_m)` with `deg P_l <= n_l`
and `P(z^s) . F(z) = P_1(z^s) f_1(z) + ... + P_m(z^s) f_m(z) = O(z^sigma)`.

`F[l][k+1]` must hold the coefficient of `z^k` in `f_l`, for `k = 0, ...,
sigma-1`. Returns

  * `P`, an `m` by `m` matrix of coefficient vectors: `P[l, lam]` holds the
    polynomial `P_{l,lam}`, so row `l` of `P` is the `l`th basis element;
  * `d`, a vector of integers with `dct P_l = d[l] + 1` (see
    [`defect`](@ref)); and
  * `pivots`, the sequence of pivot indices `pi_sigma` chosen at each step
    (`0` where no pivot was needed).

For every `delta` the corresponding solution set is

    L_delta^sigma = {a_1 P_1 + ... + a_m P_m : deg a_l <= d[l] + delta},

a `K`-vector space of dimension `sum(max(d[l] + 1 + delta, 0) for l)`.

Ties for the pivot are broken by taking the **smallest** index, which is the
uniqueness convention of the paper (§6, p. 818). If `F` has an `Integer`
element type it is promoted to `Rational{BigInt}` so the computation is
exact; the recurrence divides, so the element type must otherwise be a field.
"""
function fphps(F::AbstractVector{<:AbstractVector}, sigma::Integer,
               n::AbstractVector{<:Integer}, s::Integer)
    m = length(F)
    m >= 2 || throw(ArgumentError("need at least m = 2 series, got m=$m"))
    s >= 1 || throw(ArgumentError("s must be >= 1, got s=$s"))
    sigma >= 0 || throw(ArgumentError("sigma must be >= 0, got sigma=$sigma"))
    length(n) == m ||
        throw(ArgumentError("n must have m = $m entries, got $(length(n))"))
    for l in 1:m
        length(F[l]) >= sigma || throw(ArgumentError(
            "F[$l] must supply at least sigma = $sigma coefficients, got $(length(F[l]))"))
    end

    T = _fphps_eltype(promote_type(map(eltype, F)...))

    # P_{l,0} = e_l, d_{l,0} = n_l.
    P = Matrix{Vector{T}}(undef, m, m)
    for l in 1:m, lam in 1:m
        P[l, lam] = T[l == lam ? one(T) : zero(T)]
    end
    d = collect(Int, n)

    # R[l] carries the s-residual P_l(z^s) . F(z) truncated to degree
    # sigma-1, so that c_{l,sigma} is just R[l][sigma+1]. It is updated by
    # the same operations as P_l, which keeps each step O(sigma) instead of
    # recomputing the product from scratch.
    R = [T[k <= length(F[l]) ? F[l][k] : zero(T) for k in 1:sigma] for l in 1:m]
    pivots = Vector{Int}(undef, sigma)

    for t in 0:sigma-1
        c = [R[l][t+1] for l in 1:m]
        lambda_set = [l for l in 1:m if !iszero(c[l])]

        if isempty(lambda_set)
            pivots[t+1] = 0
            continue
        end

        # pi = the index of maximal defect in Lambda, smallest index winning
        # ties. lambda_set is ascending, so a strict > keeps the first.
        pivot = lambda_set[1]
        for l in lambda_set
            if d[l] > d[pivot]
                pivot = l
            end
        end
        pivots[t+1] = pivot

        cpi = c[pivot]
        Rp = R[pivot]
        for l in lambda_set
            l == pivot && continue
            f = c[l] / cpi
            for lam in 1:m
                _polysub!(P[l, lam], P[pivot, lam], f)
            end
            Rl = R[l]
            # Coefficients below z^t already vanish in both residuals.
            @inbounds for k in t+1:sigma
                Rl[k] -= f * Rp[k]
            end
        end

        # P_pi <- z * P_pi, hence R_pi <- z^s * R_pi.
        for lam in 1:m
            pushfirst!(P[pivot, lam], zero(T))
        end
        @inbounds for k in sigma:-1:1
            Rp[k] = k > s ? Rp[k-s] : zero(T)
        end
        d[pivot] -= 1
    end

    return P, d, pivots
end

"""
    defect(d) -> Vector{Int}

The defects `dct P_l = d[l] + 1` of the sigma-basis returned by
[`fphps`](@ref), in the sense of Beckermann and Labahn (1994), Definition
3.1: `dct P = min_l (n_l + 1 - deg P_l)`.
"""
defect(d::AbstractVector{<:Integer}) = [Int(dl) + 1 for dl in d]
