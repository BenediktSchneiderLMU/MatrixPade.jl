# MatrixPade.jl

[![CI](https://github.com/BenediktSchneiderLMU/MatrixPade.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/BenediktSchneiderLMU/MatrixPade.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/BenediktSchneiderLMU/MatrixPade.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/BenediktSchneiderLMU/MatrixPade.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

MatrixPade.jl computes Padé approximants of a **matrix-valued** power series

$$f(z) = c_0 + c_1 z + c_2 z^2 + \dots, \qquad c_k \in \mathbb{K}^{p \times q},$$

given as a plain vector of matrix coefficients. For a series of type $(M, N)$
it produces two matrix polynomials $P(z)$ (degree $\le M$) and $Q(z)$ (degree
$\le N$), with $Q$ multiplying on the right or the left, such that

$$f(z) Q(z) - P(z) = O(z^{M+N+1}), \qquad \text{or} \qquad Q(z) f(z) - P(z) = O(z^{M+N+1}).$$

This follows Beckermann & Labahn (1994), reducing the matrix problem to a
scalar simultaneous Padé (power Hermite Padé) problem, which is exposed
directly as well.

The package is pure Julia with no dependencies beyond `LinearAlgebra`, and is
generic in the coefficient number type: rational input stays rational,
symbolic input stays symbolic, and floating point is supported for the usual
numerical case.

## Installation

Not yet registered. Until then, install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/BenediktSchneiderLMU/MatrixPade.jl")
```

## Quick start

```julia
using MatrixPade

# c_0, ..., c_5 of a 2x2 matrix series (Beckermann & Labahn §5):
#
#   f(z) = [1+z^2+2z^4-z^5+...   0+...        ]
#          [0-z^5+...            1+z^2+z^4+...]
coeffs = [
    [1 0; 0 1], [0 0; 0 0], [1 0; 0 1],
    [0 0; 0 0], [2 0; 0 1], [-1 0; -1 0],
]
M, N = 2, 3

z = 1//2
f = matrix_pade(coeffs, M, N, side = :left)   # left-hand form: Q(z)^-1 P(z)
f(z)                                          # evaluate at z = 1/2, exactly
matrix_pade(coeffs, M, N, z, side = :left)    # direct evaluation at z = 1/2
```

Both evaluations give the same exact rational matrix:

```julia
2×2 Matrix{Rational{BigInt}}:
 10//7    0
 -1//21  4//3
```

The same form can be evaluated symbolically, which returns a matrix of
rational functions in `z`:

```julia
using Symbolics
@variables z
simplify.(f(z))
```

```julia
2×2 Matrix{Num}:
            (-1 - z + z^2) / (-1 - z + 2(z^2) + z^3)                0
 (z^5) / (-1 - z + 3(z^2) + 2(z^3) - 2(z^4) - (z^5))  1 / (1 - (z^2))
```

## API reference

| function | returns | notes |
|---|---|---|
| `matrix_pade_coeffs(coeffs, M, N; side, r)` | `(Pcoeffs, Qcoeffs)` | numerator first |
| `matrix_pade(coeffs, M, N; side, r)` | `MatrixPadeForm` | `side` is `:right` (default) or `:left` |
| `matrix_pade(coeffs, M, N, z; side, r)` | matrix | build and evaluate in one step |
| `matrix_pade_right(…)`, `matrix_pade_right_coeffs(…)` | as above | `side = :right` fixed |
| `matrix_pade_left(…)`, `matrix_pade_left_coeffs(…)` | as above | `side = :left` fixed |
| `MatrixPadeForm(P, Q, side)` | n/a | fields `.P`, `.Q`, `.side`; `(f)(z)` evaluates the fraction |
| `fphps(F, sigma, n, s)` | `(P, d, pivots)` | scalar σ-basis of power Hermite Padé approximants |
| `defect(d)` | `Vector{Int}` | `dct P_l = d[l] + 1` |

Coefficients are passed in mathematical order starting at $c_0$, so
`coeffs[i]` holds $c_{i-1}$. For a $p \times q$ series and type $(M, N)$:

| | defining relation | `P` | `Q` | fraction |
|---|---|---|---|---|
| right (`:right`) | `f(z) Q(z) - P(z) = O(z^{M+N+1})` | `p × r` | `q × r` | `P(z) Q(z)^-1` |
| left (`:left`) | `Q(z) f(z) - P(z) = O(z^{M+N+1})` | `r × q` | `r × p` | `Q(z)^-1 P(z)` |

`r` defaults to `q` on the right and `p` on the left, which makes `Q` square.
`fphps` is the underlying scalar solver; it solves the simultaneous-Padé
problem $P_1(z^s) f_1(z) + \dots + P_m(z^s) f_m(z) = O(z^{\sigma})$ for scalar
series and is usable on its own.

## Algorithm

Each matrix Padé form is obtained by reducing to a **scalar** power Hermite
Padé (PHPA) problem and running the FPHPS recurrence of Beckermann & Labahn
§3 (`fphps`), rather than solving the matrix problem directly:

1. The series coefficients are interleaved into a vector of `p` (or `q`)
   scalar series `F`, at stride `s = p` (right) or `s = q` (left).
2. `fphps(F, sigma, n, s)` computes a σ-basis of the solution space by
   maintaining, for each row `l`, a residual `R[l]` and a defect `d[l]`, and
   at each order reducing rows against the pivot of maximal defect. This
   keeps a single step $O(\sigma)$ and the whole run $O(m\,\sigma^2)$ instead
   of the naive cubic cost of recomputing residual products from scratch.
3. `_select_solutions` picks `r` elements out of the canonical basis of that
   solution space, preferring elements that keep the resulting `Q(0)`
   nonsingular so the denominator is invertible at the origin whenever the
   solution space allows it.
4. `P` and `Q` are read off from the selected basis elements, and evaluation
   forms `P(z) Q(z)^-1` or `Q(z)^-1 P(z)` via a pseudoinverse built from the
   normal equations (`inv` when square, `(Q'Q)\Q'`/`Q'/(Q Q')` when
   rectangular) rather than `LinearAlgebra.pinv`, which is SVD-based and
   would force floating point and destroy exactness.

### Reference

Beckermann, B. and Labahn, G. (1994). *A uniform approach for the fast
computation of matrix-type Padé approximants*. SIAM Journal on Matrix
Analysis and Applications, 15(3), 804–823.

## Limitations and assumptions

- **Square denominator is not always available.** The solution space has
  dimension at least $q(N+1) - pN$ on the right and $p(N+1) - qN$ on the
  left, so a square `Q` is guaranteed to exist only when `q >= p` on the
  right and `p >= q` on the left: a wide `f` suits `:right`, a tall `f`
  suits `:left`. Passing a smaller `r` leaves `Q` rectangular.
- **Rectangular `Q` does not fully reproduce `f`.** A rectangular `Q` has
  $QQ^{+} \ne I$, so the fraction is a valid Padé form but only reproduces
  $f(z)$ to order $M+N+1$ when `Q` is square.
- **Singular denominators are possible.** Evaluating a `MatrixPadeForm`
  throws an `ArgumentError` when $Q(z)$ is singular at the requested point.
  Beckermann & Labahn §5 exhibits a form whose $Q(z)$ is singular for
  *every* $z$ — a real feature of the theory, not an implementation gap.
- **σ-bases are not unique.** `fphps` returns *a* basis of the solution
  space; ties for the pivot are broken by the smallest index (the paper's
  own convention), but a different valid basis may differ by row
  operations.
- **Number types and exactness.** Nothing in the algorithm assumes floating
  point:

  | input `eltype` | behaviour | exact? |
  |---|---|---|
  | `Int`, `BigInt` | promoted automatically to `Rational{BigInt}` | yes |
  | `Rational{BigInt}`, `Complex{Rational{BigInt}}` | used as-is | yes |
  | `Rational{Int}` | used as-is, **overflows at modest orders** | until it overflows |
  | `Float64`, `ComplexF64` | used as-is | no |
  | `Num` (Symbolics.jl) | used as-is | yes, symbolically |

  **Use `BigInt` rationals for anything but the smallest examples.**
  Intermediate quantities grow fast; 64-bit rationals overflow sooner than
  you would guess. Plain `Int` input needs no thought, since it is promoted
  to `Rational{BigInt}` for you.
- **Out of scope.** The superfast divide-and-conquer SPHPS algorithm and the
  M-Padé generalisation to arbitrary interpolation knots of Beckermann &
  Labahn §6.

## Development

```julia
using Pkg
Pkg.test("MatrixPade")
```

or, for a faster local loop against the same tests:

```bash
julia --project=. test/runtests.jl
```

There is no build/lint step — this is a pure-Julia package with no external
dependencies beyond `LinearAlgebra` (stdlib). The suite pins the paper's
worked examples against their published values.

## License

MIT. See [LICENSE](LICENSE).
