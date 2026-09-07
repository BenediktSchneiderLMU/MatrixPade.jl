# MatrixPade.jl

[![CI](https://github.com/BenediktSchneiderLMU/MatrixPade.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/BenediktSchneiderLMU/MatrixPade.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A Julia package for evaluating matrix-Padé approximations from a vector of
matrix-valued coefficients.

> **Status:** early scaffolding. The core algorithm is not implemented yet —
> see `src/MatrixPade.jl` for the current placeholder API.

## Installation

Not yet registered. Until then, install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/BenediktSchneiderLMU/MatrixPade.jl")
```

## Usage

```julia
using MatrixPade

coeffs = [Matrix{Float64}(I, 2, 2), zeros(2, 2)]
matrix_pade_eval(coeffs, 0.5)
```

## Testing

```julia
using Pkg
Pkg.test("MatrixPade")
```

## License

MIT — see [LICENSE](LICENSE).
