# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

Run from the package root (where `Project.toml` lives).

```bash
# Full test suite (matches CI)
julia --project=. -e 'using Pkg; Pkg.test()'

# Faster local loop, with the caveat that `Symbolics` is a test-only
# dependency and so `test_symbolic.jl` errors under it
julia --project=. test/runtests.jl

# Run a single test file
julia --project=. -e 'using Test, LinearAlgebra, MatrixPade; include("test/test_matrix_pade_right.jl")'
```

There is no build/lint step — this is a pure-Julia package with no external
dependencies beyond `LinearAlgebra` (stdlib).

## Workflow

Commit regularly while developing, rather than batching everything into one
large commit at the end. A natural point to commit is once a logical unit is
done and tests pass (e.g. one algorithm piece plus its test file, one bug
fix, one doc update) — not necessarily after every single file edit.

## Architecture

**Left- and right-hand square and rectangular matrix Padé forms**, with
a genuine **matrix denominator `Q(z)`**, following Beckermann & Labahn
(1994), *"A uniform approach for the fast computation of matrix-type Padé
approximants"*, SIAM J. Matrix Anal. Appl. 15, 804–823 — their Examples 2.1
and 2.2, computed by reducing to a scalar power Hermite Padé approximation
(PHPA) problem and running the FPHPS recurrence of §3. The superfast
divide-and-conquer SPHPS variant of §6 is **out of scope**, as is the
M-Padé generalisation to arbitrary interpolation knots.

For a `p` by `q` series `A(z) = A_0 + A_1 z + ...` and a type `(M, N)`:

- right (Example 2.1): `A(z) Q(z) - P(z) = O(z^{M+N+1})`, `P` is `p` by `r`,
  `Q` is `q` by `r`, fraction `P Q^-1`.
- left (Example 2.2): `Q(z) A(z) - P(z) = O(z^{M+N+1})`, `P` is `r` by `q`,
  `Q` is `r` by `p`, fraction `Q^-1 P`.

### Source layout

- `src/fphps.jl` — the scalar solver, `fphps(F, sigma, n, s) -> (P, d,
  pivots)`. Implements the FPHPS recurrence verbatim: start from the unit
  rows `P_l = e_l` with `d_l = n_l`; at each order pick the pivot `pi` of
  maximal defect among the `l` with nonvanishing residual coefficient,
  reduce every other such row against it, then multiply row `pi` by `z` and
  drop its defect by one.
  - **Ties for `pi` are broken by the smallest index.** This is the paper's
    own "for instance" uniqueness convention (§6, p. 818) and is pinned by
    `test_fphps.jl` against the published pivot sequence of Example 4.4.
  - **The residual is carried, not recomputed.** `R[l]` holds
    `P_l(z^s) . F(z)` truncated to degree `sigma-1` and is updated by the
    same operations as `P_l` (`R_pi <- z^s R_pi` for the pivot). This is
    what keeps a step `O(sigma)` and the whole run `O(m sigma^2)`; naively
    re-forming the product would make it cubic.
  - `d` is the defect vector offset by one: `dct P_l = d[l] + 1`, exposed as
    `defect(d)`. The solution set is
    `L_delta^sigma = {sum_l a_l P_l : deg a_l <= d[l] + delta}`, of
    dimension `sum_l max(d[l] + 1 + delta, 0)`. Only `delta = 0` is used.
  - A sigma-basis is **not unique** — Definition 3.2 pins the solution
    space, not a representative — so do not assert basis matrices against
    the paper's printed ones. (Ours differs from the §5 basis on p. 815 by
    one elementary row operation.) Assert defects, orders and the derived
    `P`/`Q` instead.
- `src/matrix_pade.jl` — the two matrix Padé forms.
  - `_exactify(coeffs)` promotes `Integer`-element input matrices to
    `Rational{BigInt}` automatically so exact input stays exact
    end-to-end; non-integer input (`Float64`, `Complex`, `Rational`, ...)
    passes through unchanged. `_horner_matrix` evaluates a matrix-coefficient
    polynomial. Both live here because this is their only consumer.
  - Parameter map, Table 1 rows 2.1/2.2, verified against §5:

    | | right | left |
    |---|---|---|
    | `s` | `p` | `q` |
    | `sigma` | `p(M+N+1)` | `q(M+N+1)` |
    | `n` | `(M x p, N x q)` | `(M x q, N x p)` |
    | `F` | `F^T = (1,z,...,z^{p-1}) [I, -A(z^p)]` | `F = [I; -A(z^q)] (1,z,...,z^{q-1})^T` |
    | components | `1:p` rows of `P`, `p+1:p+q` rows of `Q`; a solution is a **column** | `1:q` cols of `P`, `q+1:q+p` cols of `Q`; a solution is a **row** |

    Concretely `f_j(z) = z^{j-1}` for `j <= s`, and the remaining entries
    interleave the coefficients of `A` at stride `s`. Only `A_0..A_{M+N}`
    are ever read.
  - `_select_solutions` picks `r` elements out of the canonical `K`-basis
    `{z^j P_l : 0 <= j <= d_l}` of `L_0^sigma`. Two passes: first greedily
    keep candidates that raise the rank of the collected `Q(0)` (so the
    denominator is nonsingular at the origin whenever the space allows it),
    then fill in canonical order. Deterministic, and it reproduces both §5
    results.
  - `r` defaults to `q` (right) / `p` (left), making `Q` square. Since
    `dim L_0^sigma >= ||n|| - sigma`, which is `q(N+1) - pN` on the right and
    `p(N+1) - qN` on the left, **a square denominator is guaranteed only
    when `q >= p` on the right and `p >= q` on the left** — a wide `A` suits
    `:right`, a tall `A` suits `:left`. The `ArgumentError` for a too-small
    solution space says so.
  - `MatrixPadeForm{T}` holds `P`, `Q` (coefficient vectors, lengths `M+1`
    and `N+1`) and `side`. `(f)(z)` evaluates `P(z) Q(z)^-1` on the right and
    `Q(z)^-1 P(z)` on the left via `_horner_matrix`.
  - `_pseudoinverse` is written by normal equations (`inv` when square,
    `(Q'Q)\Q'` when tall, `Q'/(Q Q')` when wide) rather than
    `LinearAlgebra.pinv`, which is SVD-based and would force floats and
    destroy exactness. A rectangular `Q` has `Q Q^+ != I`, so the fraction
    only reproduces `A` when `Q` is square.
  - A singular `Q(z)` raises an `ArgumentError`. This is a real case, not a
    bug: the §5 right-hand form has `Q(z)` singular for every `z` and the
    paper states no such fraction exists.
- `src/MatrixPade.jl` — module entry point; `include`s the two files above.

### Indexing convention

Series coefficients are passed as a plain `Vector` of matrices in
**mathematical order starting at `c0`**, so `coeffs[i]` (1-based Julia
index) holds `c_{i-1}`. This 0-based/1-based offset shows up throughout
`matrix_pade.jl` (e.g. `c[k+1]` for the 0-based index `k`) — keep it in mind
when touching the indexing.

### Validity constraints (throw `ArgumentError`)

`fphps`: `m >= 2`; `s >= 1`; `sigma >= 0`; `length(n) == m`; every `F[l]`
supplies at least `sigma` coefficients.

`matrix_pade_coeffs`: `side in (:right, :left)`; `M, N >= 0`;
`length(coeffs) >= M + N + 1`; all coefficient matrices the same size;
`r >= 1`; and `r <= dim L_0^sigma`.

### Tests

`test/runtests.jl` just `include`s one file per concern:

- `matrix_pade_testutils.jl` — **not** a test file; shared independent
  polynomial helpers (`_tdeg`, `_ttrim`, `_tdefect`, `_tresidual`,
  `_torder_coeff`) `include`d by the files below, written without
  package internals so the checks stay independent of the code under test.
- `exp_series_testutils.jl` — **not** a test file; builds the
  `exp(i A x + i B x^2)` series (`_texp_iab_series`, `_treadme_A`,
  `_treadme_B`) shared by the README and symbolic tests.
- `test_fphps.jl` — Beckermann & Labahn Example 4.4 (p. 813): the pivot
  sequence, the printed sigma=10 basis, the defects, and the printed
  s-residuals. Ground truth for the scalar solver alone.
- `test_matrix_pade_right.jl` — their §5 right-hand `(2,3)` form: the
  paper's exact `P` and `Q`, the defects, and the fact that evaluation must
  throw because `Q` is identically singular. Also the byproduct types
  `(0,1)` and `(1,2)` from p. 817.
- `test_matrix_pade_left.jl` — their §5 left-hand `(2,3)` form (our row
  order is the transpose-permutation of the paper's display; the form is
  not unique there), plus `Q(z) g(z) == P(z)` exactly and `O(z^{M+N+1})`
  error decay.
- `test_matrix_pade_types.jl` — exactness, shapes, the square vs.
  rectangular (pseudoinverse) paths on wide/tall series, wrapper agreement,
  and validation.
- `test_readme_example.jl` — the twelfth-order `exp(i A x + i B x^2)`
  example from the README: the series itself against its closed form, the
  order condition, exact rational evaluation, `Rational{Int}` overflow, and
  the observed `O(x^{M+N+1})` decay.
- `test_symbolic.jl` — Symbolics.jl evaluation points, checked by
  substituting a rational point back in rather than asserting a `simplify`
  normal form. Needs the test-only `Symbolics` dependency. Symbolic *series
  coefficients* are not supported: FPHPS pivots on whether a residual
  coefficient vanishes, which is not decidable for a `Num`.

When extending the algorithm, prefer encoding another worked example from
the source paper as a test the way `test_fphps.jl` and
`test_matrix_pade_right.jl` do, rather than only asserting internal
consistency.
