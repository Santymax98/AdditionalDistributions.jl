# Multivariate Distributions

AdditionalDistributions.jl provides multivariate Gaussian and Student-t distributions with support for rectangular cumulative probabilities.

A rectangular CDF is a probability of the form

```math
P(a_i \le X_i \le b_i,\quad i = 1,\ldots,d).
```

The package currently focuses on:

- `MvGaussian`, which participates directly in the
  `Distributions.AbstractMvNormal` interface;
- `MvTStudent`, a typed wrapper around `Distributions.AbstractMvTDist`;
- `cdf_result` support for both the package wrappers and native
  `Distributions.MvNormal` / `Distributions.MvTDist` objects.

```@docs
MvGaussian
MvTStudent
CDFResult
cdf_result
AdditionalDistributions.mvtcdf
```

## Basic usage

### Gaussian rectangular probabilities

```julia
using AdditionalDistributions
using LinearAlgebra
using Random

d = 5
Σ = fill(0.5, d, d)
Σ[diagind(Σ)] .= 1.0

mvnormal = MvGaussian(zeros(d), Σ)
lower = fill(-1.0, d)
upper = fill(1.0, d)

res = cdf_result(mvnormal, lower, upper;
    m = 100_000,
    rng = MersenneTwister(1234),
)

res.value
res.error
res.inform
```

`cdf(mvnormal, lower, upper)` returns only the probability estimate. Use `cdf_result` for diagnostics.

### Student-t rectangular probabilities

```julia
ν = 4.0
mvt = MvTStudent(ν, zeros(d), Σ)

res = cdf_result(mvt, lower, upper;
    m = 1_000_000,
    abseps = 1e-8,
    releps = 1e-8,
    nshifts = 16,
    rng = MersenneTwister(1234),
)
```

The Student-t CDF is usually harder than the Gaussian CDF because it uses the normal-scale-mixture representation and an additional chi-square coordinate. Small degrees of freedom, high dimension, strong dependence, and tail rectangles may require larger `m`.

## Native `Distributions.jl` interoperability

The same rectangular-probability backend is available directly for native
`Distributions.jl` multivariate distributions:

```julia
using Distributions

dn = MvNormal(zeros(d), Σ)
rn = cdf_result(dn, lower, upper;
    m = 100_000,
    rng = MersenneTwister(1234),
)

dt = MvTDist(ν, zeros(d), Σ)
rt = cdf_result(dt, lower, upper;
    m = 1_000_000,
    nshifts = 16,
    rng = MersenneTwister(1234),
)
```

This is intentionally provided through `cdf_result`, a function owned by
`AdditionalDistributions.jl`. The package does **not** add
`Distributions.cdf(::MvNormal, ...)` or `Distributions.cdf(::MvTDist, ...)`
methods, avoiding type piracy.

For package-owned wrapper types, both one-sided and rectangular scalar CDFs are
available:

```julia
cdf(mvnormal, upper)         # P(Xᵢ ≤ upperᵢ for all i)
cdf(mvnormal, lower, upper)  # P(lowerᵢ ≤ Xᵢ ≤ upperᵢ for all i)

cdf(mvt, upper)
cdf(mvt, lower, upper)
```

## API summary

```julia
cdf(dist, lower, upper; kwargs...)
```

returns only the probability estimate.

```julia
cdf_result(dist, lower, upper; kwargs...)
```

accepts `MvGaussian`, `MvTStudent`, native `MvNormal`, and native `MvTDist`
objects and returns a structured result:

```julia
res.value
res.error
res.inform
res.neval
res.algorithm
```

## Keywords

| Keyword | Meaning |
| --- | --- |
| `m` | Integration budget. |
| `abseps` | Absolute error tolerance. |
| `releps` | Relative error tolerance. |
| `rng` | Random number generator for randomized shifts. |
| `nshifts` | Number of randomized QMC shifts. Defaults: `10` for Gaussian, `16` for Student-t. |
| `batchsize` | Internal batch size. Automatic settings are usually appropriate. |
| `pivot` | Enable or disable MVSORT reordering. |
| `antithetic` | Optional antithetic reflection when available. |

Typical advanced usage:

```julia
res = cdf_result(mvt, lower, upper;
    m = 1_000_000,
    abseps = 1e-8,
    releps = 1e-8,
    nshifts = 16,
    rng = MersenneTwister(1234),
)
```

## `inform` codes

| Code | Meaning |
| ---: | --- |
| `0` | Estimated error reached the requested tolerance. |
| `1` | Estimated error is above tolerance for the current budget. |
| `2` | Invalid dimension or integration setup. |
| `3` | Matrix appears not positive semidefinite. |

`inform = 1` is common for difficult high-dimensional or heavy-tailed cases. It means the estimated integration error did not reach the requested tolerance with the current budget. Increase `m`, adjust tolerances, or compare against an external reference if the result seems suspicious.

## Algorithm summary

### Gaussian path

The Gaussian rectangular CDF uses:

1. MVSORT variable reordering;
2. Genz-style conditional transformation;
3. folded randomized Richtmyer quasi-Monte Carlo points;
4. batched evaluation to reduce per-point overhead.

Diagonal Gaussian covariance/correlation matrices are handled by an exact product shortcut.

The default number of randomized shifts is:

```julia
nshifts = 10
```

### Student-t path

The Student-t implementation uses the scale-mixture representation. If

```math
Z \sim N(0,\Sigma), \qquad W \sim \chi^2_\nu,
```

then

```math
X = \mu + \frac{Z}{\sqrt{W/\nu}}
```

has a multivariate Student-t distribution.

The algorithm uses the same conditional Gaussian core together with one radial chi-square coordinate. The Gaussian coordinates are folded; the chi-square radial coordinate is not folded.

The default settings for `MvTStudent` are more conservative:

```julia
m = max(100_000, 10_000*d)
nshifts = 16
```

Important: a diagonal Student-t scale matrix does not imply independent components. The components share a common radial scale, so diagonal Student-t rectangular probabilities are not computed as products of univariate Student-t probabilities.

## Reproducibility

Randomized QMC methods depend on random shifts. Use a fixed RNG to obtain reproducible output:

```julia
rng = MersenneTwister(1234)
res = cdf_result(mvnormal, lower, upper; rng=rng)
```

When comparing with another implementation, always report:

- `m`;
- `nshifts`;
- `abseps`;
- `releps`;
- random seed;
- full `CDFResult`;
- distribution parameters;
- lower and upper bounds.

## Reference comparisons

Gaussian reference tests include Genz-style cases and comparisons with `MvNormalCDF.jl`.

Student-t reference checks can be performed with R's `mvtnorm::pmvt` using `GenzBretz`.

See [Benchmarks](../benchmarks.md) and [Accuracy and Reproducibility](../accuracy.md) for more details.
