# Benchmarks

This page documents reproducible benchmark scripts for the multivariate rectangular CDF routines in `AdditionalDistributions.jl`.

The goal of these benchmarks is not to claim universal superiority over other implementations. Randomized quasi-Monte Carlo algorithms depend on dimension, covariance structure, integration region, integration budget, random seed, Julia version, package versions, and hardware.

Instead, these scripts provide reproducible empirical comparisons for specific configurations.

## Available scripts

Benchmark scripts live in the `benchmark/` directory.

### Basic multivariate Gaussian benchmark

```bash
julia --project=benchmark benchmark/run_mvn_basic.jl
```

This runs package-only benchmarks for `MvGaussian` rectangular probabilities.

### Comparison with `MvNormalCDF.jl`

```bash
julia --project=benchmark benchmark/run_mvn_compare_mvnormalcdf.jl
```

This benchmark compares `AdditionalDistributions.MvGaussian` against `MvNormalCDF.jl` for rectangular Gaussian probabilities under several covariance structures:

- diagonal covariance;
- equicorrelated covariance;
- AR(1)-type covariance.

The benchmark reports:

- probability estimate;
- estimated integration error;
- integration status;
- runtime;
- memory allocations.

Generated CSV files are written under:

```text
benchmark/results/
```

These generated result files should not be committed.

## Representative Gaussian results

The current `MvGaussian` implementation uses a Genz-style conditional transformation, folded randomized Richtmyer QMC points, and batched evaluation. Diagonal Gaussian rectangles are evaluated by an exact product shortcut.

In local benchmark runs, `AdditionalDistributions.jl` produced probability estimates matching `MvNormalCDF.jl` closely in the tested Gaussian cases, while often using less memory.

Representative local timings from development runs were:

| Case | AdditionalDistributions.jl | MvNormalCDF.jl | Notes |
| --- | ---: | ---: | --- |
| `d = 2`, equicorrelated `ρ = 0.2` | ≈ `0.0028s` | ≈ `0.0036s` | same probability/error order |
| `d = 5`, equicorrelated `ρ = 0.5` | ≈ `0.0155s` | ≈ `0.0179s` | same probability/error order |
| `d = 10`, equicorrelated `ρ = 0.2` | ≈ `0.0299s` | ≈ `0.0370s` | same probability/error order |
| `d = 20`, AR(1) `ρ = 0.2` | ≈ `0.0600s` | ≈ `0.0712s` | same probability/error order |

These numbers are illustrative development results, not a performance guarantee. Run the benchmark scripts on your own machine for reproducible comparisons.

## Student-t reference comparison

For `MvTStudent`, benchmark-style reference checks can be performed against R's `mvtnorm::pmvt` using the `GenzBretz` algorithm.

For example, in R:

```r
library(mvtnorm)

d <- 10
rho <- 0.5
corr <- matrix(rho, d, d)
diag(corr) <- 1

set.seed(1234)

p <- pmvt(
  lower = rep(-1, d),
  upper = rep(1, d),
  df = 4,
  corr = corr,
  algorithm = GenzBretz(
    maxpts = 1000000,
    abseps = 1e-8,
    releps = 1e-8
  )
)

print(p)
print(attributes(p))
```

A representative comparison for this case is:

```text
R mvtnorm::pmvt:
value ≈ 0.1087057
error ≈ 9.12e-6

AdditionalDistributions.jl:
value ≈ 0.1087074
error ≈ 1.46e-5
```

Both methods may report that the requested tolerance was not fully reached when strict tolerances such as `1e-8` are used.

## Reproducibility checklist

When reporting benchmark results, include:

- Julia version;
- package versions;
- operating system and CPU;
- distribution type;
- dimension;
- covariance/correlation matrix;
- integration lower and upper bounds;
- degrees of freedom for Student-t;
- `m`;
- `abseps`;
- `releps`;
- `nshifts`;
- `batchsize`;
- random seed;
- runtime;
- allocations;
- `CDFResult.value`;
- `CDFResult.error`;
- `CDFResult.inform`.

A minimal reproducible Julia example should look like:

```julia
using AdditionalDistributions
using Random, LinearAlgebra

d = 10
ρ = 0.5
Σ = fill(ρ, d, d)
Σ[diagind(Σ)] .= 1.0

lower = fill(-1.0, d)
upper = fill(1.0, d)

dist = MvGaussian(zeros(d), Σ)

res = cdf_result(dist, lower, upper;
    m = 1_000_000,
    abseps = 1e-8,
    releps = 1e-8,
    nshifts = 10,
    rng = MersenneTwister(1234),
)

res
```

For Student-t:

```julia
ν = 4.0
dist = MvTStudent(ν, zeros(d), Σ)

res = cdf_result(dist, lower, upper;
    m = 1_000_000,
    abseps = 1e-8,
    releps = 1e-8,
    nshifts = 8,
    rng = MersenneTwister(1234),
)

res
```

## Interpretation

The reported `error` is an internal randomized QMC error estimate, not a deterministic mathematical bound.

The status code `inform = 1` means that the estimated error is still above the requested tolerance for the current integration budget. It does not automatically mean that the probability estimate is unusable.

For difficult cases, increase `m`, use a fixed seed, and compare with an external reference implementation when possible.
