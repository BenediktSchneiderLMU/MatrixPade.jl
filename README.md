# MatrixPade.jl

[![CI](https://github.com/BenediktSchneiderLMU/MatrixPade.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/BenediktSchneiderLMU/MatrixPade.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/BenediktSchneiderLMU/MatrixPade.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/BenediktSchneiderLMU/MatrixPade.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Padé approximants of a **matrix-valued** power series

$$f(z) = c_0 + c_1 z + c_2 z^2 + \dots, \qquad c_k \in \mathbb{K}^{p \times q},$$

given as a plain vector of matrix coefficients. Pure Julia, no dependencies
beyond `LinearAlgebra`, and generic in the number type: rational input stays
rational, symbolic input stays symbolic.

## The two constructions

The package implements two independent algorithms from two papers. They differ
in what the denominator is, and that difference decides how many orders you get
for your coefficients.

**Matrix Padé-*type* approximants (MPTA)**, after Gu (2004). The denominator
$q(z)$ is a *scalar* polynomial of degree $n$ and the numerator $P(z)$ is
matrix-valued of degree $m$:

$$q(z) f(z) - P(z) = O(z^{m+1}).$$

**Matrix Padé forms**, after Beckermann & Labahn (1994). The denominator $Q(z)$
is a genuine *matrix* polynomial, multiplying on the right or on the left:

$$f(z) Q(z) - P(z) = O(z^{M+N+1}), \qquad Q(z) f(z) - P(z) = O(z^{M+N+1}).$$

Both consume the same number of series coefficients ($m+n+1$ and $M+N+1$), but
the MPTA matches only $m+1$ orders, while the matrix forms match all $M+N+1$.
The MPTA buys that loss back in simplicity: its "fraction" is an elementwise
divide by a scalar, so evaluating it needs no matrix inverse at all. **If you
want accuracy per coefficient, use `matrix_pade`; it is the better default.**
The MPTA is worth reaching for when you want the cheaper, structurally simpler
object, but note that its denominator does not always exist (see
[below](#limitation-the-scalar-product-can-be-isotropic)).

## Installation

Not yet registered. Until then, install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/BenediktSchneiderLMU/MatrixPade.jl")
```

## Number types

Coefficients are passed in mathematical order starting at $c_0$, so `coeffs[i]`
holds $c_{i-1}$. Nothing in either algorithm assumes floating point:

| input `eltype` | behaviour | exact? |
|---|---|---|
| `Int`, `BigInt` | promoted automatically to `Rational{BigInt}` | yes |
| `Rational{BigInt}`, `Complex{Rational{BigInt}}` | used as-is | yes |
| `Rational{Int}` | used as-is, **overflows at modest orders** | until it overflows |
| `Float64`, `ComplexF64` | used as-is | no |
| `Num` (Symbolics.jl) | used as-is | yes, symbolically |

Exactness survives evaluation too, not just the coefficient computation (here
`coeffs` is the integer series of the next section):

```julia
julia> matrix_pade_left(coeffs, 2, 3, 1//2)    # exact rational in, exact rational out
2×2 Matrix{Rational{BigInt}}:
  10//7   0//1
 -1//21   4//3
```

The linear algebra is written to keep this property. In particular the
pseudoinverse used for rectangular denominators is formed from the normal
equations rather than `LinearAlgebra.pinv`, which is SVD-based and would force
floats.

**Use `BigInt` rationals for anything but the smallest examples.** The
intermediate determinants and cofactors grow fast, and 64-bit rationals overflow
sooner than you would guess: for the twelfth-order example below,
`Complex{Rational{Int}}` already dies at $(4,4)$:

```
ERROR: OverflowError: -23434417 * 3093361140689 overflowed for type Int64
```

whereas `Complex{Rational{BigInt}}` runs cleanly through $(6,6)$. Plain `Int`
input needs no thought, since it is promoted to `Rational{BigInt}` for you.

## Matrix Padé forms

For a $p \times q$ series $f(z)$, a matrix Padé form of type $(M, N)$ is a pair
of matrix polynomials with $\deg P \le M$, $\deg Q \le N$ satisfying an order
condition. There are two, according to which side the denominator multiplies on:

| | defining relation | `P` | `Q` | fraction |
|---|---|---|---|---|
| right (`:right`) | `f(z) Q(z) - P(z) = O(z^{M+N+1})` | `p × r` | `q × r` | `P(z) Q(z)^-1` |
| left (`:left`) | `Q(z) f(z) - P(z) = O(z^{M+N+1})` | `r × q` | `r × p` | `Q(z)^-1 P(z)` |

Only $c_0, \dots, c_{M+N}$ are used. `r` defaults to `q` on the right and `p` on
the left, making `Q` square. Since the solution space has dimension at least
$q(N+1) - pN$ on the right and $p(N+1) - qN$ on the left, **a square denominator
is guaranteed to exist only when `q >= p` on the right and `p >= q` on the
left**: a wide `f` suits `:right`, a tall `f` suits `:left`. Passing a smaller
`r` leaves `Q` rectangular, and evaluation then uses the corresponding one-sided
Moore–Penrose pseudoinverse; note that a rectangular `Q` has $QQ^{+} \ne I$, so
the fraction only reproduces $f(z)$ to order $M+N+1$ when `Q` is square.

```julia
using MatrixPade

# The example of Beckermann & Labahn §5, as coefficients c_0, ..., c_5 of
#
#   f(z) = [1+z^2+2z^4-z^5+...   0+...        ]
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
x = 1//2
g(x)                                  # evaluate at z = 1/2, exactly

# ...or the general entry point, and immediate evaluation:
matrix_pade(coeffs, 2, 3; side=:left)
matrix_pade(coeffs, 2, 3, x; side=:left)
```

`f.P` and `f.Q` are the numerator and denominator coefficient lists, so
`P(z) = f.P[1] + f.P[2] z + ...`. Evaluation throws an `ArgumentError` if `Q(z)`
is singular at the requested point. For the example above that is exactly what
happens on the right, where `Q(z)` is singular for *every* `z` and no right
matrix Padé fraction of type $(2,3)$ exists: a real feature of the theory, not
a defect of the implementation.

The underlying scalar solver is exposed as `fphps(F, sigma, n, s)`, which
returns a σ-basis of power Hermite Padé approximants together with the defect
vector (see `defect`) and the pivot sequence. It solves the simultaneous-Padé
problem $P_1(z^s) f_1(z) + \dots + P_m(z^s) f_m(z) = O(z^{\sigma})$ for scalar
series, and is useful on its own.

## Matrix Padé-type approximants

```julia
using MatrixPade

# f(z) = exp(A*z), A = [0 1; 0 -2], as a matrix power series c0 + c1 z + ...
coeffs = [
    [1//1 0//1; 0//1 1//1],
    [0//1 1//1; 0//1 -2//1],
    [0//1 -1//1; 0//1 2//1],
    [0//1 2//3; 0//1 -4//3],
    [0//1 -1//3; 0//1 2//3],
    [0//1 2//15; 0//1 -4//15],
]

# Coefficient lists directly (note: denominator first):
qcoeffs, Pcoeffs = mpta_coeffs(coeffs, 3, 2)   # q(z), P(z) for the (3/2) MPTA

# ...or as a callable approximant:
pa = mpta(coeffs, 3, 2)
x = 1//2
pa(1//2)                                        # evaluate at z = 1/2

# ...or evaluate directly:
mpta(coeffs, 3, 2, x)
```

This is Example 4.4 of Gu (2004), and the package reproduces the published
$q$ and $P$ exactly: $\det H = 50/27$ and
$q(z) = \tfrac{50}{27} + \tfrac{40}{27} z + \tfrac{10}{27} z^2$.

### The approximant is a rational function, so ask for it symbolically

Because everything is generic in the number type, you can evaluate at a
*symbolic* point and read off the closed form:

```julia
using Symbolics
@variables z

julia> simplify.(mpta(coeffs, 3, 2)(z))
2×2 Matrix{Num}:
    1           (z*((50//27) + z*(-(10//27) + (10//81)*z))) / ((50//27) + z*((40//27) + (10//27)*z))
 0//1  ((50//27) + z*(-(20//9) + z*((10//9) - (20//81)*z))) / ((50//27) + z*((40//27) + (10//27)*z))
```

Clearing denominators, this is

$$(3/2)_f(z) = \begin{pmatrix} 1 & \dfrac{z(15 - 3z + z^2)}{3(5 + 4z + z^2)} \\[2ex] 0 & \dfrac{15 - 18z + 9z^2 - 2z^3}{15 + 12z + 3z^2} \end{pmatrix},$$

which is the paper's result. Note that the $(1,1)$ and $(2,1)$ entries come out
as exactly $1$ and $0$, not merely correct to some order, because
$\exp(Az)$ has those entries identically for this $A$.

### Limitation: the scalar product can be isotropic

`matdot(A, B) = sum(A .* B)` is the *bilinear* product Gu's construction is
built on. It deliberately does **not** conjugate, unlike the usual Hermitian
inner product. On a complex series that means it can vanish on nonzero
coefficients, and if enough of them are isotropic the Hankel matrix is singular
and no Padé-type denominator exists. The package detects this and throws rather
than returning a `NaN`-valued approximant:

```
ERROR: ArgumentError: the Padé-type denominator is identically zero: the Hankel
matrix of matdot products is singular for (m, n) = (3, 3). …
```

The example in the next section is exactly such a case, and it is not a corner
case: *every* diagonal $(k,k)$ approximant of it degenerates. When this
happens, use `matrix_pade`, whose matrix denominator does not depend on
`matdot` at all.

## Worked example: $\exp(iAx + iBx^2)$ to twelfth order

Take

$$A = \begin{pmatrix} 1 & i \\ -i & -1 \end{pmatrix}, \qquad B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \qquad M(x) = i\left(A x + B x^2\right),$$

and approximate $f(x) = e^{M(x)}$.

$A$ and $B$ anticommute, so $M(x)$ at one value of $x$ does not commute with
$M(x)$ at another: the eigenbasis of $M(x)$ *depends on $x$*. The standard
shortcut for $e^{Ax}$ (diagonalise once, build scalar Padé approximants of the
eigenvalues, rotate back) does not apply here, so the approximation has to be
done on the matrix series itself.

There is still an exact reference to check against, which is what makes the
example useful. Since $A^2 = 2I$, $B^2 = I$ and $AB = -BA$,

$$M(x)^2 = -x^2(2 + x^2)\, I \quad \Longrightarrow \quad e^{M(x)} = \cos\theta \, I + \frac{\sin\theta}{\theta} M(x), \qquad \theta = x\sqrt{2 + x^2},$$

which is the exact closed form used as the reference below.

### Building the series

```julia
using MatrixPade, LinearAlgebra

T = Complex{Rational{BigInt}}
A = T[1 im; -im -1]
B = T[0 1; 1 0]
N = 12

# Coefficients of exp(M(x)) by truncated exponentiation of M(x) = i(A x + B x^2).
function exp_series(A, B, N)
    T = eltype(A); d = size(A, 1)
    zed() = zeros(T, d, d); eye() = Matrix{T}(I, d, d)
    M = [zed(), im * A, im * B]
    function tmul(a, b)                       # product truncated at degree N
        r = [zed() for _ in 0:N]
        for i in eachindex(a), j in eachindex(b)
            k = (i - 1) + (j - 1)
            k <= N && (r[k+1] += a[i] * b[j])
        end
        r
    end
    S = [zed() for _ in 0:N]; S[1] = eye()
    Mk = vcat([eye()], [zed() for _ in 1:N]); fact = big(1)
    for k in 1:N
        Mk = tmul(Mk, M); fact *= k
        for j in 0:N; S[j+1] += Mk[j+1] / fact; end
    end
    S
end

coeffs = exp_series(A, B, N)     # c_0, ..., c_12, exactly
```

The $(1,1)$ entries of $c_0, \dots, c_{12}$ come out as

$$1,\; i,\; -1,\; -\tfrac{i}{3},\; -\tfrac13,\; -\tfrac{2i}{15},\; \tfrac{7}{45},\; \tfrac{2i}{63},\; \tfrac{8}{315},\; \tfrac{17i}{2835},\; -\tfrac{107}{14175},\; -\tfrac{172i}{155925},\; -\tfrac{109}{133650},$$

which agree exactly with the closed form $e^{M(x)}$ above, and likewise for the
other three entries.

### Exact vs. approximated

Thirteen coefficients are exactly what the $(6,6)$ form needs
($M + N + 1 = 13$), so the whole ladder from $(1,1)$ to $(6,6)$ is available:

```julia
cf = [ComplexF64.(c) for c in coeffs]
Af, Bf = ComplexF64.(A), ComplexF64.(B)
exact(x) = exp(im * Af * x + im * Bf * x^2)

x = 0.5
for k in 1:6
    approx = matrix_pade_right(cf, k, k, x)
    println(k, "  ", norm(approx - exact(x)))
end
```

`norm(approx - exact)` at each order:

| $(M,N)$ | $x = 0.5$ | $x = 1.0$ |
|---|---|---|
| $(1,1)$ | 1.505e-01 | 9.497e-01 |
| $(2,2)$ | 5.692e-03 | 1.735e-01 |
| $(3,3)$ | 5.711e-04 | 5.910e-02 |
| $(4,4)$ | 5.816e-06 | 3.184e-03 |
| $(5,5)$ | 6.200e-07 | 1.059e-03 |
| $(6,6)$ | 2.581e-09 | 2.333e-05 |

`matrix_pade_left` agrees with `matrix_pade_right` to the digits shown. The
$O(x^{M+N+1})$ rate is visible by halving the argument, $\mathrm{err}(x) /
\mathrm{err}(x/2)$ at $x = 0.1$:

| $(M,N)$ | $(1,1)$ | $(2,2)$ | $(3,3)$ | $(4,4)$ |
|---|---|---|---|---|
| observed ratio | 7.97 | 32.04 | 127.83 | 510.37 |
| $2^{M+N+1}$ | 8 | 32 | 128 | 512 |


### The same thing symbolically and exactly

Evaluate the $(2,2)$ form at a symbolic point and it hands back the rational
function itself. With $D(x) = 18 + 3x^2 + 5x^4$:

$$\left[\tfrac{2}{2}\right]_f(x) = \frac{1}{D(x)}\begin{pmatrix} 18 - 15x^2 - 4x^4 + i(18x - 3x^3) & -18x + 3x^3 + i(18x^2 - 3x^4) \\[1ex] 18x - 3x^3 + i(18x^2 - 3x^4) & 18 - 15x^2 - 4x^4 - i(18x - 3x^3) \end{pmatrix}$$

Substituting $x = 1/2$ gives $D = 305/16$ and a $(1,1)$ entry of
$(224 + 138i)/305$, which is exactly what evaluating the approximant on
rationals returns, with no floating point anywhere:

```julia
julia> matrix_pade_right(coeffs, 2, 2, 1//2)
2×2 Matrix{Complex{Rational{BigInt}}}:
 224//305+138//305*im  -138//305+69//305*im
 138//305+69//305*im    224//305-138//305*im
```

The $(1,1)$ entry along the whole ladder at $x = 1/2$ shows how quickly the
exact values outgrow 64-bit rationals:

| $(M,N)$ | `matrix_pade_right(coeffs, k, k, 1//2)[1,1]` |
|---|---|
| $(1,1)$ | `4//5 + 2//5*im` |
| $(2,2)$ | `224//305 + 138//305*im` |
| $(3,3)$ | `46864//64025 + 29082//64025*im` |
| $(4,4)$ | `268089248//366396473 + 166499290//366396473*im` |
| $(5,5)$ | `21949205668979//29997991820821 + 13631856216360//29997991820821*im` |
| $(6,6)$ | `6608375858586849397//9031674702833386045 + 4104226355542786536//9031674702833386045*im` |

## API reference

| function | returns | notes |
|---|---|---|
| `mpta_coeffs(coeffs, m, n)` | `(qcoeffs, Pcoeffs)` | **denominator first** |
| `mpta(coeffs, m, n)` | `PadeApproximant` | callable; `pa(z)` gives `P(z) ./ q(z)` |
| `mpta(coeffs, m, n, z)` | matrix | build and evaluate in one step |
| `PadeApproximant(q, P)` | n/a | fields `.q` (scalars), `.P` (matrices) |
| `matdot(A, B)` | scalar | `sum(A .* B)`; bilinear, **does not conjugate** |
| `matrix_pade_coeffs(coeffs, M, N; side, r)` | `(Pcoeffs, Qcoeffs)` | **numerator first** |
| `matrix_pade(coeffs, M, N; side, r)` | `MatrixPadeForm` | `side` is `:right` (default) or `:left` |
| `matrix_pade(coeffs, M, N, z; side, r)` | matrix | build and evaluate |
| `matrix_pade_right(…)`, `matrix_pade_right_coeffs(…)` | as above | `side = :right` fixed |
| `matrix_pade_left(…)`, `matrix_pade_left_coeffs(…)` | as above | `side = :left` fixed |
| `MatrixPadeForm(P, Q, side)` | n/a | fields `.P`, `.Q`, `.side` |
| `fphps(F, sigma, n, s)` | `(P, d, pivots)` | scalar σ-basis of power Hermite Padé approximants |
| `defect(d)` | `Vector{Int}` | `dct P_l = d[l] + 1` |

**Watch the return order.** `mpta_coeffs` returns `(q, P)`: denominator first,
while `matrix_pade_coeffs` returns `(P, Q)`: numerator first. Each follows the
convention of its own source paper.

## Limitations

- **MPTA on complex series.** `matdot` is bilinear and non-conjugating, so it
  can be isotropic; the denominator then does not exist and `mpta_coeffs`
  throws. See above.
- **Singular denominators.** `MatrixPadeForm` evaluation throws an
  `ArgumentError` when $Q(z)$ is singular at the requested point. Beckermann &
  Labahn §5 exhibits a form whose $Q(z)$ is singular for *every* $z$, so this is
  a genuine possibility, not an implementation gap.
- **Rectangular denominators.** With `r` smaller than the default, $Q$ is
  rectangular, $QQ^{+} \ne I$, and the fraction no longer reproduces $f(z)$ to
  order $M+N+1$. It is still a valid Padé *form*.
- **σ-bases are not unique.** `fphps` returns *a* basis of the solution space;
  ties for the pivot are broken by the smallest index, which is the paper's own
  convention, but a different valid basis may differ by row operations.
- **Out of scope.** The directional MPTA and the system-reduction application
  of Gu §§5–6; the superfast SPHPS algorithm and the M-Padé generalisation to
  arbitrary interpolation knots of Beckermann & Labahn §6.

## Testing

```julia
using Pkg
Pkg.test("MatrixPade")
```

The suite pins both papers' worked examples against their published values, and
covers the twelfth-order example and the symbolic paths documented above.

## Citing

If you use this package, please cite the underlying papers:

```bibtex
@article{Gu2004,
  author  = {Gu, Chuanqing},
  title   = {Matrix {P}ad{\'e}-type approximant and directional matrix
             {P}ad{\'e} approximant in the inner product space},
  journal = {Journal of Computational and Applied Mathematics},
  volume  = {164--165},
  pages   = {365--385},
  year    = {2004}
}

@article{BeckermannLabahn1994,
  author  = {Beckermann, Bernhard and Labahn, George},
  title   = {A uniform approach for the fast computation of matrix-type
             {P}ad{\'e} approximants},
  journal = {SIAM Journal on Matrix Analysis and Applications},
  volume  = {15},
  number  = {3},
  pages   = {804--823},
  year    = {1994}
}
```

## License

MIT. See [LICENSE](LICENSE).
