# AdditionalDistributions.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Santymax98.github.io/AdditionalDistributions.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Santymax98.github.io/AdditionalDistributions.jl/dev/)
[![Build Status](https://github.com/Santymax98/AdditionalDistributions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Santymax98/AdditionalDistributions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Santymax98/AdditionalDistributions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Santymax98/AdditionalDistributions.jl)

**AdditionalDistributions.jl** extends the [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl) ecosystem with additional continuous, discrete, and multivariate probability distributions.

The package follows the standard `Distributions.jl` interface whenever possible: `pdf`, `logpdf`, `cdf`, `quantile`, `rand`, `mean`, `var`, `params`, `minimum`, `maximum`, and related statistical methods.

It also provides native Julia routines for rectangular cumulative probabilities of multivariate Gaussian and Student-t distributions.

---

## Installation

AdditionalDistributions.jl is registered in the Julia General registry.

```julia
using Pkg
Pkg.add("AdditionalDistributions")
```

Then load it with:

```julia
using AdditionalDistributions
```

To install the development version directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/Santymax98/AdditionalDistributions.jl")
```

---

## Highlights

- Additional continuous and discrete univariate distributions.
- Zero-inflated discrete models such as `ZIP`, `ZIB`, and `ZINB`.
- Heavy-tailed distributions such as `Burr`, `Lomax`, and `Zeta`.
- Rectangular CDF evaluation for `MvGaussian`.
- Rectangular CDF evaluation for `MvTStudent`.
- Structured multivariate CDF diagnostics through `CDFResult` and `cdf_result`.
- Reproducible randomized QMC integration when a seeded RNG is provided.
- Reference tests and benchmark scripts for multivariate CDF evaluation.

---

## Basic usage

```julia
using AdditionalDistributions
using Distributions

x = Lomax(2.0, 3.0)
pdf(x, 1.5)
cdf(x, 1.5)

z = ZIP(2.0, 0.3)
pdf(z, 0)
cdf(z, 4)
```

Additional distributions are documented in the [stable documentation](https://santymax98.github.io/AdditionalDistributions.jl/stable/) and in the [distribution index](https://santymax98.github.io/AdditionalDistributions.jl/stable/Distributions/).

---

## Multivariate rectangular probabilities

The multivariate API computes probabilities of rectangles

```math
P(a_i \le X_i \le b_i,\quad i = 1,\ldots,d).
```

```julia
using AdditionalDistributions
using LinearAlgebra
using Random

d = 10
Σ = fill(0.5, d, d)
Σ[diagind(Σ)] .= 1.0

lower = fill(-1.0, d)
upper = fill( 1.0, d)

mvnormal = MvGaussian(zeros(d), Σ)

res = cdf_result(mvnormal, lower, upper;
    m = 100_000,
    rng = MersenneTwister(1234),
)

res.value      # estimated probability
res.error      # estimated integration error
res.inform     # integration status code
res.neval      # integration budget
res.algorithm  # algorithm identifier
```

`cdf(dist, lower, upper)` returns only the probability estimate. Use `cdf_result` when you need diagnostics.

For Student-t:

```julia
ν = 4.0
mvt = MvTStudent(ν, zeros(d), Σ)

res_t = cdf_result(mvt, lower, upper;
    m = 1_000_000,
    abseps = 1e-8,
    releps = 1e-8,
    nshifts = 16,
    rng = MersenneTwister(1234),
)
```

The multivariate Student-t implementation uses the normal-scale-mixture representation. Even with a diagonal scale matrix, its components are not treated as independent univariate Student-t random variables, because they share a common radial scale.

---

## Multivariate CDF controls

| Keyword | Meaning |
| --- | --- |
| `m` | Maximum integration budget. Larger values usually reduce QMC noise. |
| `abseps` | Absolute error tolerance. |
| `releps` | Relative error tolerance. |
| `rng` | Random number generator used for randomized shifts. Fix it for reproducibility. |
| `nshifts` | Number of randomized QMC shifts. Defaults: `12` for `MvGaussian`, `16` for `MvTStudent`. |
| `batchsize` | Internal batch size. `nothing` or an automatic setting lets the implementation choose a value. |
| `pivot` | Whether to use the MVSORT variable reordering step. |
| `antithetic` | Optional antithetic reflection when available. |

`inform` codes:

| Code | Meaning |
| ---: | --- |
| `0` | Estimated error is within tolerance. |
| `1` | Estimated error is above tolerance for the current budget. Increase `m` or relax tolerances. |
| `2` | Invalid dimension or integration setup. |
| `3` | Matrix appears not positive semidefinite during preparation. |

`inform = 1` does **not** mean that the probability estimate is invalid. It means that the internal error estimator did not meet the requested tolerance with the current integration budget.

---

## Benchmarks

Benchmark scripts are available under `benchmark/`.

```bash
julia --project=benchmark benchmark/run_mvn_basic.jl
julia --project=benchmark benchmark/run_mvn_compare_mvnormalcdf.jl
```

The comparison benchmark checks `MvGaussian` against `MvNormalCDF.jl` for diagonal, equicorrelated, and AR(1)-type covariance structures.

Generated CSV files are written to `benchmark/results/` and are not tracked by Git.

Benchmark results depend on the random seed, dimension, covariance structure, integration bounds, integration budget, Julia version, package versions, and machine. They should be interpreted as reproducible empirical comparisons for specific configurations, not as universal performance guarantees.

See the [benchmarks documentation](https://santymax98.github.io/AdditionalDistributions.jl/stable/benchmarks/) for details.

---

## Running tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Useful local filters are documented in `test/runtests.jl`, for example:

```julia
@run_package_tests (filter = ti -> :multivariate in ti.tags)
@run_package_tests (filter = ti -> :continuous in ti.tags)
@run_package_tests (filter = ti -> :discrete in ti.tags)
```

---

## Reporting numerical issues

Please include:

- distribution type (`MvGaussian` or `MvTStudent`),
- dimension,
- `μ`, `Σ`, lower bounds and upper bounds,
- degrees of freedom `ν` for Student-t,
- `m`, `abseps`, `releps`, `nshifts`, `batchsize`, `pivot`,
- RNG type and seed,
- full `CDFResult(value, error, inform, neval, algorithm)`,
- Julia version and platform,
- external reference result, if available.

---

## Contributing

Contributions are welcome. Please see [`CONTRIBUTING.md`](CONTRIBUTING.md) for development notes, test conventions, benchmark guidelines, and how to report reproducible numerical cases.

Do not hesitate to star this repository to show support.

---

## Citation

If you use AdditionalDistributions.jl in your research, please cite it as:

> S. Jiménez (2025). *AdditionalDistributions.jl — Advanced and Extended Probability Distributions in Julia*. Available at: https://github.com/Santymax98/AdditionalDistributions.jl

BibTeX:

```bibtex
@misc{Jimenez2025AdditionalDistributions,
    author = {Santiago Jimenez},
    title = {AdditionalDistributions.jl --- Advanced and Extended Probability Distributions in Julia},
    year = {2025},
    url = {https://github.com/Santymax98/AdditionalDistributions.jl},
    note = {Julia package}
}
```

---

## License

MIT License.
