# MatrixPade.jl

[![CI](https://github.com/BenediktSchneiderLMU/MatrixPade.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/BenediktSchneiderLMU/MatrixPade.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A Julia package for computing Padé approximants of a matrix-valued power
series, given as a vector of matrix coefficients.

Currently implemented:

- The `(m/n)` **matrix Padé-type approximant** (MPTA) with a **scalar**
  denominator `q(z)` and a **matrix-valued** numerator `P(z)`, satisfying
  `q(z) f(z) - P(z) = O(z^{m+1})`. Follows Gu (2004), *"Matrix Padé-type
  approximant and directional matrix Padé approximant in the inner product
  space"*, J. Comput. Appl. Math. 164–165, 365–385.
- The **right-** and **left-hand square and rectangular matrix Padé forms**,
  with a **matrix** denominator `Q(z)`. Follows Beckermann & Labahn (1994),
  *"A uniform approach for the fast computation of matrix-type Padé
  approximants"*, SIAM J. Matrix Anal. Appl. 15, 804–823 (Examples 2.1 and
  2.2, computed with the FPHPS algorithm of §3).

In both cases, integer-valued input is promoted to exact `Rational{BigInt}`
automatically, so exact input gives exact coefficients and exact evaluation;
`Float64`/`Complex` input is left as-is.

## Installation

Not yet registered. Until then, install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/BenediktSchneiderLMU/MatrixPade.jl")
```

## Matrix Padé forms

For a `p` by `q` matrix power series `A(z) = A_0 + A_1 z + ...`, a matrix
Padé form of type `(M, N)` is a pair of matrix polynomials with
`deg P <= M`, `deg Q <= N` satisfying an order condition. There are two,
according to which side the denominator multiplies on:

| | defining relation | `P` | `Q` | fraction |
|---|---|---|---|---|
| right (`:right`) | `A(z) Q(z) - P(z) = O(z^{M+N+1})` | `p` by `r` | `q` by `r` | `P(z) Q(z)^-1` |
| left (`:left`) | `Q(z) A(z) - P(z) = O(z^{M+N+1})` | `r` by `q` | `r` by `p` | `Q(z)^-1 P(z)` |

Only `A_0, ..., A_{M+N}` are used. `r` defaults to `q` on the right and `p`
on the left, making `Q` square; a square denominator is guaranteed to exist
only when `q >= p` on the right and `p >= q` on the left, so a wide `A`
suits `:right` and a tall `A` suits `:left`. Passing a smaller `r` leaves
`Q` rectangular, and evaluation then uses the corresponding one-sided
Moore–Penrose pseudoinverse (computed by normal equations, so exact input
stays exact).

```julia
using MatrixPade

# The example of Beckermann & Labahn §5, as coefficients A_0, ..., A_5 of
#
#   A(z) = [1+z^2+2z^4-z^5+...   0+...        ]
#          [0-z^5+...            1+z^2+z^4+...]
coeffs = [
    [1 0; 0 1], [0 0; 0 0], [1 0; 0 1],
    [0 0; 0 0], [2 0; 0 1], [-1 0; -1 0],
]

# Coefficient lists directly:
Pcoeffs, Qcoeffs = matrix_pade_right_coeffs(coeffs, 2, 3)  # P(z), Q(z)

# ...or as a callable form evaluating P(z) * Q(z)^-1:
f = matrix_pade_right(coeffs, 2, 3)

# ...the left-hand form instead, evaluating Q(z)^-1 * P(z):
g = matrix_pade_left(coeffs, 2, 3)
g(1//2)                                  # evaluate at z = 1/2, exactly

# ...or the general entry point, and immediate evaluation:
matrix_pade(coeffs, 2, 3; side=:left)
matrix_pade(coeffs, 2, 3, 1//2; side=:left)
```

`f.P` and `f.Q` are the numerator and denominator coefficient lists, so
`P(z) = f.P[1] + f.P[2] z + ...`. Evaluation throws an `ArgumentError` if
`Q(z)` is singular at the requested point — for the example above that is
exactly what happens on the right, where `Q(z)` is singular for every `z`
and no right matrix Padé *fraction* of type `(2, 3)` exists.

The underlying scalar solver is exposed as `fphps(F, sigma, n, s)`, which
returns a sigma-basis of power Hermite Padé approximants together with the
defect vector and the pivot sequence.

## Matrix Padé-type approximants

```julia
using MatrixPade

# f(z) = exp(A*z), A = [0 1; 0 -2], as a matrix power series c0 + c1 z + ...
coeffs = [
    [1 0; 0 1],
    [0 1; 0 -2],
    [0 -1; 0 2],
    [0 2//3; 0 -4//3],
    [0 -1//3; 0 2//3],
    [0 2//15; 0 -4//15],
]

# Coefficient lists directly:
qcoeffs, Pcoeffs = mpta_coeffs(coeffs, 3, 2)   # q(z), P(z) for the (3/2) MPTA

# ...or as a callable approximant:
pa = mpta(coeffs, 3, 2)
pa(1//2)                                        # evaluate at z = 1/2

# ...or evaluate directly:
mpta(coeffs, 3, 2, 1//2)
```

## Testing

```julia
using Pkg
Pkg.test("MatrixPade")
```

## License

MIT — see [LICENSE](LICENSE).
