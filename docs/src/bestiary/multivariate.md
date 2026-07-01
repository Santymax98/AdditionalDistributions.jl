# Multivariate Distributions

AdditionalDistributions.jl provides wrappers around selected `Distributions.jl` multivariate laws with rectangular CDF support.

A rectangular CDF is a probability of the form

```math
P(a_i \le X_i \le b_i,\quad i = 1,\ldots,d).
```

The package currently focuses on:

- `MvGaussian`, based on `Distributions.MvNormal`;
- `MvTStudent`, based on `Distributions.MvTDist`.

```@docs
MvGaussian
MvTStudent
CDFResult
cdf_result
AdditionalDistributions.mvtcdf
```

## Gaussian rectangular probabilities

```julia
using AdditionalDistributions
using LinearAlgebra
using Random

d = 5
Σ = fill(0.5, d, d)
Σ[diagind(Σ)] .= 1.0

mvnormal = MvGaussian(zeros(d), Σ)
a = fill(-1.0, d)
b = fill(1.0, d)

res = cdf_result(mvnormal, a, b;
    m = 100_000,
    rng = MersenneTwister(1234),
)

res.value
res.error
res.inform
```

`cdf(mvnormal, a, b)` returns only the probability. Use `cdf_result` for diagnostics.

## Student's t rectangular probabilities

```julia
ν = 4.0
mvt = MvTStudent(ν, zeros(d), Σ)

res = cdf_result(mvt, a, b;
    m = 1_000_000,
    abseps = 1e-8,
    releps = 1e-8,
    nshifts = 16,
    rng = MersenneTwister(1234),
)
```

The Student's t CDF is usually harder than the Gaussian CDF because it uses the normal-scale-mixture representation and an additional chi-square coordinate. Small degrees of freedom, high dimension, and high correlations may require larger `m`.

## Algorithm summary

The Gaussian path uses:

1. MVSORT variable reordering;
2. Genz-style conditional transformation;
3. folded randomized Richtmyer quasi-Monte Carlo points;
4. batch evaluation to reduce per-point overhead.

The Student's t path uses the same conditional Gaussian core together with the scale-mixture representation. The Gaussian coordinates are folded; the chi-square radial coordinate is not folded.

## Keywords

| Keyword | Meaning |
| --- | --- |
| `m` | Integration budget. |
| `abseps` | Absolute error tolerance. |
| `releps` | Relative error tolerance. |
| `rng` | Random number generator for randomized shifts. |
| `nshifts` | Number of randomized QMC shifts. Defaults: `12` for Gaussian, `16` for Student's t. |
| `batchsize` | Internal batch size. `0` selects an automatic setting. |
| `pivot` | Enable/disable MVSORT reordering. |
| `antithetic` | Optional antithetic reflection in Gaussian coordinates. |

## `inform` codes

| Code | Meaning |
| --- | --- |
| `0` | Estimated error reached the requested tolerance. |
| `1` | Estimated error is above tolerance for the current budget. |
| `2` | Invalid dimension. |
| `3` | Matrix appears not positive semidefinite. |

`inform = 1` is common for difficult high-dimensional or heavy-tailed cases. Increase `m`, adjust tolerances, or report the full reproducible case if the result seems suspicious.

## Reproducibility

Randomized QMC methods depend on random shifts. Use a fixed RNG to obtain reproducible output:

```julia
rng = MersenneTwister(1234)
res = cdf_result(mvnormal, a, b; rng=rng)
```

When comparing with another implementation, always report `m`, `nshifts`, `abseps`, `releps`, the seed, and the full `CDFResult`.
