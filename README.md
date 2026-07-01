# AdditionalDistributions.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Santymax98.github.io/AdditionalDistributions.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Santymax98.github.io/AdditionalDistributions.jl/dev/)
[![Build Status](https://github.com/Santymax98/AdditionalDistributions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Santymax98/AdditionalDistributions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Santymax98/AdditionalDistributions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Santymax98/AdditionalDistributions.jl)

**AdditionalDistributions.jl** extends the [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl) ecosystem with additional discrete, continuous, and multivariate probability distributions.

The package follows the standard `Distributions.jl` interface whenever possible: `pdf`, `logpdf`, `cdf`, `quantile`, `rand`, `mean`, `var`, `params`, and related statistical methods.

## Highlights

- Additional **continuous** and **discrete** univariate distributions, including zero-inflated and heavy-tailed families.
- `MvGaussian`, a wrapper around `Distributions.MvNormal` with rectangular CDF evaluation.
- `MvTStudent`, a wrapper around `Distributions.MvTDist` with rectangular CDF evaluation.
- Structured multivariate CDF output through `cdf_result` and `CDFResult`.
- Pure-Julia randomized quasi-Monte Carlo integration for multivariate rectangular probabilities.
- Reference tests against known formulas, Genz-style multivariate normal cases, and selected `mvtnorm::pmvt` values from R.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Santymax98/AdditionalDistributions.jl")
```

## Basic use

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

## Multivariate rectangular probabilities

The multivariate API computes probabilities of rectangles

```math
P(a_i \le X_i \le b_i,\; i=1,\ldots,d).
```

```julia
using AdditionalDistributions
using LinearAlgebra
using Random

d = 10
Σ = fill(0.5, d, d)
Σ[diagind(Σ)] .= 1.0

a = fill(-1.0, d)
b = fill( 1.0, d)

mvnormal = MvGaussian(zeros(d), Σ)
p = cdf(mvnormal, a, b; m=100_000, rng=MersenneTwister(1234))

res = cdf_result(mvnormal, a, b;
    m = 100_000,
    rng = MersenneTwister(1234),
)

res.value      # estimated probability
res.error      # estimated absolute integration error
res.inform     # convergence/status code
res.algorithm  # algorithm identifier
```

For Student's t:

```julia
ν = 4.0
mvt = MvTStudent(ν, zeros(d), Σ)
res_t = cdf_result(mvt, a, b;
    m = 1_000_000,
    abseps = 1e-8,
    releps = 1e-8,
    nshifts = 16,
    rng = MersenneTwister(1234),
)
```

## Multivariate CDF controls

| Keyword | Meaning |
| --- | --- |
| `m` | Maximum integration budget. Larger values usually reduce QMC noise. |
| `abseps` | Absolute error tolerance. |
| `releps` | Relative error tolerance. |
| `rng` | Random number generator used for randomized shifts. Fix it for reproducibility. |
| `nshifts` | Number of randomized QMC shifts. Defaults: `12` for `MvGaussian`, `16` for `MvTStudent`. |
| `batchsize` | Internal batch size. `0` selects an automatic memory/speed tradeoff. |
| `pivot` | Whether to use the MVSORT variable reordering step. |
| `full=true` | Legacy option returning `(value, error, inform)` from `cdf`. |

`inform` codes:

| Code | Meaning |
| --- | --- |
| `0` | Estimated error is within tolerance. |
| `1` | Estimated error is above tolerance for the current budget. Increase `m` or relax tolerances. |
| `2` | Invalid dimension. |
| `3` | Matrix appears not positive semidefinite during preparation. |

`inform = 1` does **not** mean that the probability is unusable. It means that the internal error estimator did not meet the requested tolerance with the current budget.

## Numerical method

The multivariate normal CDF uses a Genz-style conditional transformation with MVSORT reordering and folded randomized Richtmyer quasi-Monte Carlo evaluation in batches. Independent Gaussian rectangles are evaluated exactly by factorization.

The multivariate Student's t CDF uses the scale-mixture representation of the t distribution, applying the same conditional Gaussian transformation with an additional radial chi-square coordinate. The chi-square coordinate is not folded.

## Benchmarks and references

The benchmark script

```bash
julia --project=benchmark benchmark/run_mvn_compare_mvnormalcdf.jl
```

compares `MvGaussian` against `MvNormalCDF.jl` for diagonal, equicorrelated, and AR(1) covariance structures. Generated CSV files are intentionally not tracked.

For Student's t, selected values are checked against R's `mvtnorm::pmvt` using `GenzBretz`. Because randomized QMC estimates depend on the seed, always report the seed and integration settings when comparing results.

## Running tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Useful local filters are available in `test/runtests.jl`, for example:

```julia
@run_package_tests (filter = ti -> :multivariate in ti.tags)
@run_package_tests (filter = ti -> :continuous in ti.tags)
@run_package_tests (filter = ti -> :discrete in ti.tags)
```

## Reporting numerical issues

Please include:

- distribution type (`MvGaussian` or `MvTStudent`),
- dimension,
- `μ`, `Σ`, lower bounds and upper bounds,
- degrees of freedom `ν` for Student's t,
- `m`, `abseps`, `releps`, `nshifts`, `batchsize`, `pivot`,
- RNG and seed,
- `CDFResult(value, error, inform, neval, algorithm)`,
- Julia version and platform,
- external reference result, if available.

## Contributing

Contributions are welcome. Please see [`CONTRIBUTING.md`](CONTRIBUTING.md) for development notes, test conventions, benchmark guidelines, and how to report reproducible numerical cases.

## License

MIT License.
