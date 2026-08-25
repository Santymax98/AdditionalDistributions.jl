# Accuracy and Reproducibility

`AdditionalDistributions.jl` uses randomized quasi-Monte Carlo methods for multivariate rectangular probabilities. The returned `error` is an internal estimate based on randomized shifts, not a rigorous deterministic bound.

This page explains how to interpret numerical output, how to reproduce results, and what information to include when reporting numerical issues.

## `CDFResult`

For multivariate rectangular CDFs, prefer `cdf_result` when you need diagnostic information:

```julia
res = cdf_result(dist, lower, upper;
    m = 100_000,
    abseps = 1e-6,
    releps = 1e-6,
    rng = MersenneTwister(1234),
)
```

The returned object stores:

| Field | Meaning |
| --- | --- |
| `value` | Estimated rectangular probability. |
| `error` | Estimated randomized QMC integration error. |
| `inform` | Integration status code. |
| `neval` | Requested integration budget. The actual prime-rounded lattice evaluation count may be slightly smaller. |
| `algorithm` | Internal algorithm identifier. |

## `inform` codes

| Code | Meaning |
| ---: | --- |
| `0` | Estimated error is within the requested tolerance. |
| `1` | Estimated error is above tolerance for the current budget. |
| `2` | Invalid dimension or integration setup. |
| `3` | Matrix appears not positive semidefinite during preparation. |

`inform = 1` is not an exception and does not automatically invalidate the returned probability. It means that the estimated error did not reach `abseps`/`releps` with the current integration budget.

## Choosing `m`

The integration budget `m` controls the approximate number of function evaluations. Larger values are typically required for:

- larger dimension;
- stronger correlation;
- tail rectangles;
- narrow rectangles;
- Student-t distributions with small degrees of freedom.

For quick exploratory computations, moderate values of `m` may be enough. For publication-quality numerical comparisons, use a fixed seed and increase `m` until the estimate is stable for the precision you need.

## `nshifts`

`nshifts` controls the number of randomized QMC shifts used to estimate uncertainty.

Defaults:

- `MvGaussian`: `nshifts = 10`;
- `MvTStudent`: `nshifts = 8`.

For fixed `m`, increasing `nshifts` reduces the number of points per shift. It does not always reduce the reported error. Benchmark before changing defaults.

## Reproducible results

Always set the random number generator when you need reproducibility:

```julia
using Random

rng = MersenneTwister(1234)

res = cdf_result(dist, lower, upper;
    m = 1_000_000,
    abseps = 1e-8,
    releps = 1e-8,
    rng = rng,
)
```

Small differences across seeds are expected because the method is randomized.

## Reference comparisons

For structured Gaussian validation cases, deterministic references are used whenever available. Equicorrelated rectangles are reduced to one-dimensional quadrature, while AR(1) cases are evaluated through a deterministic Gaussian Markov recursion.

For Student-t validation cases, the normal-scale-mixture representation is combined with the same deterministic Gaussian references and a one-dimensional radial integral.

`MvNormalCDF.jl`, SciPy, and R's `mvtnorm` are useful independent implementation comparisons, but their randomized estimates are not treated as ground truth.

For difficult Student-t cases, especially with small degrees of freedom, high dimension, or strong correlation, larger integration budgets may still be required.

## Recommended reporting format

When reporting a numerical issue, include a complete reproducible example:

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
nshifts = 10

# For Student-t:
# dist = MvTStudent(ν, zeros(d), Σ)
# nshifts = 8

res = cdf_result(dist, lower, upper;
    m = 1_000_000,
    abseps = 1e-8,
    releps = 1e-8,
    nshifts = nshifts,
    rng = MersenneTwister(1234),
)

res
```

Then report:

- `res.value`;
- `res.error`;
- `res.inform`;
- `res.neval`;
- `res.algorithm`;
- distribution parameters and bounds;
- Julia version and platform;
- package versions;
- external reference values, if available.
