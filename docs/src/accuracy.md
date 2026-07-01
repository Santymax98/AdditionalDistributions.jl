# Accuracy and Reproducibility

This package uses randomized quasi-Monte Carlo methods for multivariate rectangular probabilities. The returned `error` is an internal estimate based on randomized shifts, not a rigorous deterministic bound.

## Recommended reporting format

When reporting a numerical issue, include:

```julia
using AdditionalDistributions
using Random

res = cdf_result(d, a, b;
    m = 100_000,
    abseps = 1e-6,
    releps = 1e-6,
    nshifts = 12,
    batchsize = 0,
    rng = MersenneTwister(1234),
)
```

Then report:

- `res.value`,
- `res.error`,
- `res.inform`,
- `res.neval`,
- `res.algorithm`,
- the distribution parameters and bounds,
- Julia version and platform,
- external reference values if available.

## Choosing `m`

The integration budget `m` controls the approximate number of function evaluations. Larger values are typically required for:

- larger dimension,
- stronger correlation,
- tail rectangles,
- Student's t distributions with small degrees of freedom.

## `nshifts`

`nshifts` controls the number of randomized QMC shifts used to estimate uncertainty.

Defaults:

- `MvGaussian`: `nshifts = 12`;
- `MvTStudent`: `nshifts = 16`.

For fixed `m`, increasing `nshifts` reduces the number of points per shift. It does not always reduce the reported error. Benchmark before changing defaults.

## Reference comparisons

Gaussian cases are tested against Genz-style reference values adapted from the `MvNormalCDF.jl` test suite.

Student's t cases are tested against selected values generated with R's `mvtnorm::pmvt` and the `GenzBretz` algorithm.
