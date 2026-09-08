# Left- and right-hand square and rectangular matrix Padé forms:
# Beckermann and Labahn (1994), Examples 2.1 and 2.2. Both are reduced to a
# scalar PHPA problem via the parameter map of Table 1 (rows 2.1 and 2.2)
# and solved with the FPHPS recurrence of §3.

_coeff(poly::AbstractVector{T}, j::Integer) where {T} =
    j + 1 <= length(poly) ? poly[j+1] : zero(T)

# F^T(z) = (1, z, ..., z^{p-1}) . [I, -A(z^p)], i.e. f_j(z) = z^{j-1} for
# j <= p and f_{p+j}(z) = -sum_i z^{i-1} A_{ij}(z^p) for j = 1, ..., q.
function _build_F_right(c, M, N, p, q, sigma, ::Type{T}) where {T}
    F = [zeros(T, sigma) for _ in 1:(p+q)]
    for j in 1:p
        j <= sigma && (F[j][j] = one(T))
    end
    for j in 1:q, i in 1:p, k in 0:(M+N)
        idx = (i - 1) + p * k + 1
        idx <= sigma || continue
        F[p+j][idx] -= c[k+1][i, j]
    end
    return F
end

# F(z) = [I; -A(z^q)] . (1, z, ..., z^{q-1})^T, i.e. f_j(z) = z^{j-1} for
# j <= q and f_{q+i}(z) = -sum_j A_{ij}(z^q) z^{j-1} for i = 1, ..., p.
function _build_F_left(c, M, N, p, q, sigma, ::Type{T}) where {T}
    F = [zeros(T, sigma) for _ in 1:(p+q)]
    for j in 1:q
        j <= sigma && (F[j][j] = one(T))
    end
    for i in 1:p, j in 1:q, k in 0:(M+N)
        idx = (j - 1) + q * k + 1
        idx <= sigma || continue
        F[q+i][idx] -= c[k+1][i, j]
    end
    return F
end

# Reduce v against the already-collected rows; append it and return true if
# it is linearly independent of them.
function _add_independent!(rows, pivots, v::AbstractVector{T}) where {T}
    w = collect(T, v)
    for (row, pv) in zip(rows, pivots)
        iszero(w[pv]) || (w -= (w[pv] / row[pv]) * row)
    end
    pv = findfirst(!iszero, w)
    pv === nothing && return false
    push!(rows, w)
    push!(pivots, pv)
    return true
end

# The K-basis of L_0^sigma is {z^j P_l : 0 <= j <= d[l]}. Pick `r` of those
# elements, preferring ones that raise the rank of the collected Q(0) so the
# denominator is nonsingular at the origin whenever the space allows it,
# then filling up in canonical order.
function _select_solutions(P, d, r, qrange, ::Type{T}) where {T}
    m = size(P, 1)
    cands = [(l, j) for l in 1:m for j in 0:d[l]]
    length(cands) >= r || throw(ArgumentError(
        "cannot extract r = $r solutions: the solution space has dimension " *
        "$(length(cands)); try a smaller r or larger M, N"))

    chosen = Int[]
    rows = Vector{Vector{T}}()
    pivots = Int[]
    for (idx, (l, j)) in enumerate(cands)
        length(chosen) < r || break
        # The constant term of z^j P_{l,lam} vanishes unless j == 0.
        v = [j == 0 ? _coeff(P[l, lam], 0) : zero(T) for lam in qrange]
        _add_independent!(rows, pivots, v) && push!(chosen, idx)
    end
    for idx in eachindex(cands)
        length(chosen) < r || break
        idx in chosen || push!(chosen, idx)
    end
    sort!(chosen)

    return [(l = cands[idx][1], shift = cands[idx][2]) for idx in chosen]
end

# Coefficient of z^j in the lam-th component of the selected solution
# z^shift * P_l.
_solution_coeff(P, sol, lam, j, ::Type{T}) where {T} =
    j < sol.shift ? zero(T) : _coeff(P[sol.l, lam], j - sol.shift)

"""
    matrix_pade_coeffs(coeffs, M, N; side=:right, r=nothing) -> (Pcoeffs, Qcoeffs)

Compute a matrix Padé form of type `(M, N)` for the matrix power series
`A(z) = coeffs[1] + coeffs[2] z + ...` with `A` of size `p` by `q`,
following Beckermann and Labahn (1994).

For `side = :right` (their Example 2.1) the result satisfies

    A(z) Q(z) - P(z) = z^{M+N+1} R(z),   P is p by r,  Q is q by r,

and for `side = :left` (their Example 2.2)

    Q(z) A(z) - P(z) = z^{M+N+1} R(z),   P is r by q,  Q is r by p,

with `deg P <= M` and `deg Q <= N`. Returns `Pcoeffs = [P_0, ..., P_M]` and
`Qcoeffs = [Q_0, ..., Q_N]`, so `P(z) = P_0 + P_1 z + ... + P_M z^M`.

`r` defaults to `q` on the right and `p` on the left, which makes `Q` square;
passing a smaller `r` gives a rectangular denominator, which is still a valid
Padé *form* but no longer determines `A` (see [`MatrixPadeForm`](@ref)).
A square denominator is guaranteed to exist only when `q >= p` on the right
and `p >= q` on the left, so a wide `A` suits `:right` and a tall `A` suits
`:left`; otherwise an `ArgumentError` reports the dimension actually
available.

Only the first `M+N+1` series coefficients are used. If the coefficient
matrices have an `Integer` element type they are promoted to
`Rational{BigInt}` so the computation is exact.

See [`matrix_pade`](@ref) for the callable form, and [`fphps`](@ref) for the
underlying scalar solver.
"""
function matrix_pade_coeffs(coeffs::AbstractVector{<:AbstractMatrix}, M::Integer, N::Integer;
                            side::Symbol=:right, r::Union{Nothing,Integer}=nothing)
    side in (:right, :left) ||
        throw(ArgumentError("side must be :right or :left, got :$side"))
    M >= 0 || throw(ArgumentError("M must be >= 0, got M=$M"))
    N >= 0 || throw(ArgumentError("N must be >= 0, got N=$N"))
    length(coeffs) >= M + N + 1 || throw(ArgumentError(
        "need at least M+N+1 = $(M + N + 1) coefficients, got $(length(coeffs))"))
    p, q = size(coeffs[1])
    all(size(A) == (p, q) for A in coeffs) || throw(ArgumentError(
        "all coefficient matrices must have the same size, got $(unique(size.(coeffs)))"))

    c = _exactify(coeffs)
    T = eltype(c[1])

    if side === :right
        s, sigma = p, p * (M + N + 1)
        n = vcat(fill(Int(M), p), fill(Int(N), q))
        F = _build_F_right(c, M, N, p, q, sigma, T)
        qrange = (p+1):(p+q)
        rr = r === nothing ? q : Int(r)
    else
        s, sigma = q, q * (M + N + 1)
        n = vcat(fill(Int(M), q), fill(Int(N), p))
        F = _build_F_left(c, M, N, p, q, sigma, T)
        qrange = (q+1):(q+p)
        rr = r === nothing ? p : Int(r)
    end
    rr >= 1 || throw(ArgumentError("r must be >= 1, got r=$rr"))

    B, d, _ = fphps(F, sigma, n, s)

    # dim L_0^sigma >= ||n|| - sigma, which is q(N+1) - pN on the right and
    # p(N+1) - qN on the left. So a square denominator is guaranteed to exist
    # only when q >= p on the right, and when p >= q on the left.
    dim = sum(max(dl + 1, 0) for dl in d)
    dim >= rr || throw(ArgumentError(
        "cannot extract r = $rr solutions from a solution space of dimension " *
        "$dim for the $side-hand form of a $p by $q series of type ($M, $N); " *
        "pass a smaller r, increase N, or use side = " *
        ":$(side === :right ? "left" : "right") (the $side-hand form needs " *
        "$(side === :right ? "q >= p" : "p >= q") for a square denominator)"))

    sols = _select_solutions(B, d, rr, qrange, T)

    if side === :right
        # Each solution is a column: components 1:p give P, p+1:p+q give Q.
        Pcoeffs = [T[_solution_coeff(B, sols[k], i, j, T) for i in 1:p, k in 1:rr]
                   for j in 0:M]
        Qcoeffs = [T[_solution_coeff(B, sols[k], p + i, j, T) for i in 1:q, k in 1:rr]
                   for j in 0:N]
    else
        # Each solution is a row: components 1:q give P, q+1:q+p give Q.
        Pcoeffs = [T[_solution_coeff(B, sols[k], i, j, T) for k in 1:rr, i in 1:q]
                   for j in 0:M]
        Qcoeffs = [T[_solution_coeff(B, sols[k], q + i, j, T) for k in 1:rr, i in 1:p]
                   for j in 0:N]
    end

    return Pcoeffs, Qcoeffs
end

"""
    MatrixPadeForm(P, Q, side)

A matrix Padé form with matrix-valued numerator coefficients `P` and
denominator coefficients `Q`, as produced by
[`matrix_pade_coeffs`](@ref). `side` is `:right` or `:left`.

Calling `f(z)` evaluates the rational function that approximates `A(z)`:
`P(z) * Q(z)^-1` for `:right` and `Q(z)^-1 * P(z)` for `:left`. When `Q(z)`
is rectangular the corresponding one-sided Moore-Penrose pseudoinverse is
used; note that the fraction reproduces `A(z)` to order `M+N+1` only when
`Q(z)` is square (the default), since a rectangular `Q` has `Q Q^+ != I`.

Evaluation throws an `ArgumentError` if `Q(z)` is singular at the requested
point. That is a genuine possibility rather than a defect: §5 of the paper
exhibits a right-hand form whose `Q(z)` is singular for every `z`, so no
matrix Padé fraction of that type exists.
"""
struct MatrixPadeForm{T<:AbstractMatrix}
    P::Vector{T}
    Q::Vector{T}
    side::Symbol
end

# Moore-Penrose pseudoinverse for the shapes this package produces, written
# so that exact (Rational) input stays exact -- LinearAlgebra.pinv is
# SVD-based and would force floating point.
function _pseudoinverse(Q::AbstractMatrix, z)
    a, b = size(Q)
    try
        if a == b
            return inv(Q)
        elseif a > b
            return (Q' * Q) \ Matrix(Q')   # full column rank
        else
            return Matrix(Q') / (Q * Q')   # full row rank
        end
    catch e
        e isa SingularException || e isa LinearAlgebra.LAPACKException || rethrow()
        throw(ArgumentError(
            "the denominator Q(z) is singular at z = $z, so the matrix Padé " *
            "fraction does not exist there"))
    end
end

function (f::MatrixPadeForm)(z)
    Pz = _horner_matrix(f.P, z)
    Qz = _horner_matrix(f.Q, z)
    Qinv = _pseudoinverse(Qz, z)
    return f.side === :right ? Pz * Qinv : Qinv * Pz
end

"""
    matrix_pade(coeffs, M, N; side=:right, r=nothing) -> MatrixPadeForm
    matrix_pade(coeffs, M, N, z; side=:right, r=nothing) -> matrix

Build the type `(M, N)` matrix Padé form of `A(z) = coeffs[1] + coeffs[2] z +
...` (see [`matrix_pade_coeffs`](@ref)) as a [`MatrixPadeForm`](@ref). If `z`
is given, evaluate it immediately and return the resulting matrix instead.

[`matrix_pade_right`](@ref) and [`matrix_pade_left`](@ref) are convenience
wrappers that fix `side`.
"""
matrix_pade(coeffs::AbstractVector{<:AbstractMatrix}, M::Integer, N::Integer;
            side::Symbol=:right, r::Union{Nothing,Integer}=nothing) =
    MatrixPadeForm(matrix_pade_coeffs(coeffs, M, N; side=side, r=r)..., side)

matrix_pade(coeffs::AbstractVector{<:AbstractMatrix}, M::Integer, N::Integer, z;
            side::Symbol=:right, r::Union{Nothing,Integer}=nothing) =
    matrix_pade(coeffs, M, N; side=side, r=r)(z)

"""
    matrix_pade_right_coeffs(coeffs, M, N; r=nothing) -> (Pcoeffs, Qcoeffs)
    matrix_pade_right(coeffs, M, N; r=nothing) -> MatrixPadeForm
    matrix_pade_right(coeffs, M, N, z; r=nothing) -> matrix

Right-hand matrix Padé form of type `(M, N)`, satisfying
`A(z) Q(z) - P(z) = z^{M+N+1} R(z)` and approximating `A(z)` by
`P(z) Q(z)^-1`. Shorthand for [`matrix_pade`](@ref) with `side = :right`.
"""
matrix_pade_right_coeffs(coeffs::AbstractVector{<:AbstractMatrix}, M::Integer, N::Integer;
                         r::Union{Nothing,Integer}=nothing) =
    matrix_pade_coeffs(coeffs, M, N; side=:right, r=r)

matrix_pade_right(coeffs::AbstractVector{<:AbstractMatrix}, M::Integer, N::Integer;
                  r::Union{Nothing,Integer}=nothing) =
    matrix_pade(coeffs, M, N; side=:right, r=r)

matrix_pade_right(coeffs::AbstractVector{<:AbstractMatrix}, M::Integer, N::Integer, z;
                  r::Union{Nothing,Integer}=nothing) =
    matrix_pade(coeffs, M, N, z; side=:right, r=r)

"""
    matrix_pade_left_coeffs(coeffs, M, N; r=nothing) -> (Pcoeffs, Qcoeffs)
    matrix_pade_left(coeffs, M, N; r=nothing) -> MatrixPadeForm
    matrix_pade_left(coeffs, M, N, z; r=nothing) -> matrix

Left-hand matrix Padé form of type `(M, N)`, satisfying
`Q(z) A(z) - P(z) = z^{M+N+1} R(z)` and approximating `A(z)` by
`Q(z)^-1 P(z)`. Shorthand for [`matrix_pade`](@ref) with `side = :left`.
"""
matrix_pade_left_coeffs(coeffs::AbstractVector{<:AbstractMatrix}, M::Integer, N::Integer;
                        r::Union{Nothing,Integer}=nothing) =
    matrix_pade_coeffs(coeffs, M, N; side=:left, r=r)

matrix_pade_left(coeffs::AbstractVector{<:AbstractMatrix}, M::Integer, N::Integer;
                 r::Union{Nothing,Integer}=nothing) =
    matrix_pade(coeffs, M, N; side=:left, r=r)

matrix_pade_left(coeffs::AbstractVector{<:AbstractMatrix}, M::Integer, N::Integer, z;
                 r::Union{Nothing,Integer}=nothing) =
    matrix_pade(coeffs, M, N, z; side=:left, r=r)
